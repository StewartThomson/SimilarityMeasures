using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Statistics;

namespace SimilarityMeasures{
    public static class SimilarityMeasures{

        private static double MinOf3(double a, double b, double c){
            return Math.Min(a, Math.Min(b, c));
        }
        private static bool TrajCheck(Matrix<double> trajectory1, Matrix<double> trajectory2){
            if(trajectory1.ColumnCount != trajectory2.ColumnCount){
                Console.WriteLine("Trajectory dimensions do not match");
                return false;
            }

            return true;
        }

        public static Vector<double> AveTranslate(Matrix<double> trajectory1, Matrix<double> trajectory2){
            if(!TrajCheck(trajectory1, trajectory2)){
                return null;
            }
            
            int dimensions = trajectory1.ColumnCount;
            int length1 = trajectory1.RowCount;
            int length2 = trajectory2.RowCount;

            Vector<double> translation = Vector<double>.Build.Dense(dimensions, 0.0);

            if(dimensions == 0){
                Console.WriteLine("Dimension is 0");
                return translation;
            }

            if(length1 == 0 || length2 == 0){
                Console.WriteLine("At least one trajectory contains 0 points");
                return translation;
            }

            double[][] columns1 = trajectory1.ToColumnArrays();
            double[][] columns2 = trajectory2.ToColumnArrays();
            for(int i = 0; i < dimensions; i++){
                double mean1 = ArrayStatistics.Mean(columns1[i]);
                double mean2 = ArrayStatistics.Mean(columns2[i]);;
                double newTranslation = mean1 - mean2;
                translation[i] = newTranslation;
            }

            return translation;
        }

        public static bool DistanceCheck(Vector<double> point1, Vector<double> point2, double distance){
            int dimensions = point1.Count;
            bool check = true;

            for(int i = 0; i < dimensions; i++){
                double newDist = Math.Abs(point1[i] - point2[i]);

                if(newDist > distance){
                    check = false;
                    break;
                }
            }
            return check;
        }

        public static double DistanceSq(Vector<double> point1, Vector<double> point2){
            int dimensions = point1.Count;

            double dist = 0;

            for(int i = 0; i < dimensions; i++){
                dist += Math.Pow(point1[i] - point2[i], 2);
            }

            return dist;
        }

        public static double DynamicTimeWarping(Matrix<double> trajectory1, Matrix<double> trajectory2, int pointSpacing = -1){
            if(!TrajCheck(trajectory1, trajectory2)){
                return -1;
            }

            int dimensions = trajectory1.ColumnCount;
            int length1 = trajectory1.RowCount;
            int length2 = trajectory2.RowCount;

            if(length1 == 0 || length2 == 0){
                Console.WriteLine("At least one trajectory contains 0 points");
                return -1;
            }
            if(dimensions == 0){
                Console.WriteLine("Dimension is 0");
                return -1;
            }

            //Default point spacing
            if(pointSpacing < 0){
                pointSpacing = Math.Max(length1, length2);
            }

            //Length1 rows, length2 columns, populated with -1
            Matrix<double> warpPaths = Matrix<double>.Build.Dense(length1, length2, -1);
            double dist = DistanceSq(trajectory1.Row(1), trajectory2.Row(1));
            warpPaths[0,0] = Math.Sqrt(dist);

            //Initializing matrices
            if(length1 > 1 & pointSpacing > 0){
                for(int i = 1; i < Math.Min(length1, pointSpacing + 1); i++){
                    dist = DistanceSq(trajectory1.Row(i), trajectory2.Row(1));
                    warpPaths[i, 0] = Math.Sqrt(dist) + warpPaths[i - 1, 0];
                }
            }

            if(length2 > 1 & pointSpacing > 0){
                for(int i = 1; i < Math.Min(length2, pointSpacing + 1); i++){
                    dist = DistanceSq(trajectory1.Row(1), trajectory2.Row(i));
                    warpPaths[0, i] = Math.Sqrt(dist) + warpPaths[0, i - 1];
                }
            }

            //Set up rest of warp path matrix
            if(length1 > 1 & length2 > 1 & pointSpacing >= 0){
                for(int point1 = 1; point1 < length1; point1++){
                    for(int point2 = 1; point2 < length2; point2++){
                        int pointDifference = point1 - point2;

                        //When within point distance
                        if(Math.Abs(pointDifference) <= pointSpacing){
                            dist = DistanceSq(trajectory1.Row(point1), trajectory2.Row(point2));
                            double path = -1;
                            if(pointSpacing == 0){
                                path = warpPaths[point1 - 1, point2 - 1];
                            }else if(pointDifference == pointSpacing){
                                path = Math.Min(warpPaths[point1 - 1, point2 -1 ], warpPaths[point1 - 1, point2]);
                            }else if(-pointDifference == pointSpacing){
                                path = Math.Min(warpPaths[point1 - 1, point2 - 1], warpPaths[point1, point2 - 1]);
                            }else{
                                path = Math.Min(warpPaths[point1 - 1, point2 - 1], Math.Min(warpPaths[point1, point2 - 1], warpPaths[point1 -1, point2]));
                            }

                            warpPaths[point1, point2] = path + Math.Sqrt(dist);
                        }
                    }
                }
            }

            return warpPaths[length1 - 1, length2 - 1];
        }

        public static int EditDist(Matrix<double> trajectory1, Matrix<double> trajectory2, double pointDistance = 20){
            if(!TrajCheck(trajectory1, trajectory2)){
                return -1;
            }

            int dimensions = trajectory1.ColumnCount;
            int length1 = trajectory1.RowCount;
            int length2 = trajectory2.RowCount;

            if(length1 == 0){
                Console.WriteLine("Trajectory 1 has no points");
                return length2;
            }
            if(length2 == 0){
                Console.WriteLine("Trajectory 2 has no points");
                return length1;
            }

            if(dimensions == 0){
                Console.WriteLine("Dimension is 0");
                return -1;
            }

            Matrix<double> editPaths = Matrix<double>.Build.Dense(length1 + 1, length2 + 1, -1);
            for(int i = 0; i < length1 + 1; i++){
                 editPaths[i, 0] = i;
            }

            for(int i = 1; i < length2 + 1; i++){
                editPaths[0, i] = i;
            }

            for(int point1 = 1; point1 < length1 + 1; point1++){
                for(int point2 = 1; point2 < length2 + 1; point2++){
                    int diagonal = 2;

                    if(DistanceCheck(trajectory1.Row(point1 - 1), trajectory2.Row(point2 - 1), pointDistance)){
                        diagonal = 1;
                    }

                    double pathValue = MinOf3(editPaths[point1 - 1, point2] + 1, editPaths[point1, point2 - 1] + 1, editPaths[point1 - 1, point2 - 1] + diagonal);

                    editPaths[point1, point2] = pathValue;
                }
            }
            Console.WriteLine(editPaths);
            return (int)editPaths[length1, length2];
        }
    }    
}