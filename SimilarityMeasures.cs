using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Statistics;
using System.Collections.Generic;

namespace SimilarityMeasures{
    public static class SimilarityMeasures{

        private static double MinOf3(double a, double b, double c){
            return Math.Min(a, Math.Min(b, c));
        }
        /*Checks that the trajectories have the same dimension
        Args: 
            trajectory1: An m x n matrix containing trajectory1. Here m is the number 
            of points and n is the dimension of the points.
            trajectory2: traj2: A k x n matrix containing trajectory2. Here k is the number 
            of points and n is the dimension of the points.
            The two trajectories do not need to have the same number of points
        Returns:
            True if the dimensions are the same, false otherwise
         */
        private static bool TrajCheck(Matrix<double> trajectory1, Matrix<double> trajectory2){
            if(trajectory1.ColumnCount != trajectory2.ColumnCount){
                Console.WriteLine("Trajectory dimensions do not match");
                return false;
            }

            return true;
        }

        /*Calculates the translation vector for trajectory 2 using the average of the two trajectories
        Args: 
            trajectory1: An m x n matrix containing trajectory1. Here m is the number 
            of points and n is the dimension of the points.
            trajectory2: traj2: A k x n matrix containing trajectory2. Here k is the number 
            of points and n is the dimension of the points.
            The two trajectories do not need to have the same number of points
        Returns:
            A vector of length n is returned containing the translation in each dimension. If there is a problem an empty vector is returned.
         */
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

        /*Checks whether two points lie within some distance in every dimension
        Args:
            point1: An n dimensional vector representing point1
            point2: An n dimensional vector representing point2
            dist: Number representing the maximum distance in each dimension allowed for points to be considered equivalent
        Returns:
            True if the points are within the distance, false if not
         */
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

        /*Calculates the square of the distance between two points with the same dimensions
        Args:
            point1: An n dimensional vector representing point1
            point2: An n dimensional vector representing point2
        Returns:
            A double representing the square of the distance between the two points
         */
        public static double DistanceSq(Vector<double> point1, Vector<double> point2){
            int dimensions = point1.Count;

            double dist = 0;

            for(int i = 0; i < dimensions; i++){
                dist += Math.Pow(point1[i] - point2[i], 2);
            }

            return dist;
        }

        /*Calculates the Dynamic Time Warping value between two trajectories
        Args: 
            trajectory1: An m x n matrix containing trajectory1. Here m is the number 
                of points and n is the dimension of the points.
            trajectory2: traj2: A k x n matrix containing trajectory2. Here k is the number 
                of points and n is the dimension of the points.
                The two trajectories do not need to have the same number of points
            pointSpacing: An integer value of the maximum index difference between trajectory1 and
                trajectory2 allowed in the calculation. A negative values ets the spacing to unlimited
        Returns:
            A number representing the smallest warp path is returned. If a problem occurs, -1 is returned
        */
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
                            //When no point spacing is allowed
                            if(pointSpacing == 0){
                                path = warpPaths[point1 - 1, point2 - 1];
                            }else if(pointDifference == pointSpacing){
                                //The furthest distance forward point calculation
                                path = Math.Min(warpPaths[point1 - 1, point2 -1 ], warpPaths[point1 - 1, point2]);
                            }else if(-pointDifference == pointSpacing){
                                //The furthest distance backwards point calculation
                                path = Math.Min(warpPaths[point1 - 1, point2 - 1], warpPaths[point1, point2 - 1]);
                            }else{
                                //All other points
                                path = Math.Min(warpPaths[point1 - 1, point2 - 1], Math.Min(warpPaths[point1, point2 - 1], warpPaths[point1 -1, point2]));
                            }
                            //Storing the new best path for these points
                            warpPaths[point1, point2] = path + Math.Sqrt(dist);
                        }
                    }
                }
            }
            //Returning the best warp path distance
            return warpPaths[length1 - 1, length2 - 1];
        }

        /*Calculates the edit distance between two trajectories.
        This takes two trajectories and the maximum distance between points
        Args: 
            trajectory1: An m x n matrix containing trajectory1. Here m is the number 
                of points and n is the dimension of the points.
            trajectory2: traj2: A k x n matrix containing trajectory2. Here k is the number 
                of points and n is the dimension of the points.
                The two trajectories do not need to have the same number of points
            pointDistance: A number representing the maximum distance in each dimension is allowed
                for points to be considered equivalent
        Returns:
            An integer representing the minimum number of edits required to be returned. If a problem occurs
            -1 is returned.
         */
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
                    int diagonal = 1;
                    //Setting the diagonal increment depending on whether the current points are within range
                    if(DistanceCheck(trajectory1.Row(point1 - 1), trajectory2.Row(point2 - 1), pointDistance)){
                        diagonal = 0;
                    }
                    //Setting the path to the current two points to the minimum value
                    double pathValue = MinOf3(editPaths[point1 - 1, point2] + 1, editPaths[point1, point2 - 1] + 1, editPaths[point1 - 1, point2 - 1] + diagonal);

                    editPaths[point1, point2] = pathValue;
                }
            }
            //Returning the final minimum path as the edit distance
            return (int)editPaths[length1, length2];
        }

        /*Calculates the longest common subsequence for two given trajectories
        Args: 
            trajectory1: An m x n matrix containing trajectory1. Here m is the number 
                of points and n is the dimension of the points.
            trajectory2: traj2: A k x n matrix containing trajectory2. Here k is the number 
                of points and n is the dimension of the points.
                The two trajectories do not need to have the same number of points
            pointSpacing: An integer value of the maximum index difference between trajectory1 and
                trajectory2 allowed in the calculation. A negative values ets the spacing to unlimited
            pointDistance: A number representing the maximum distance in each dimension is allowed
                for points to be considered equivalent
            errorMargin: A number used to scale the accuracy and speed of the calculation
        Returns:
            A vector containing the LCSS value (in its first index) and the translations are retuned.
            The first value of the vector is the LCSS values. If a problem occurs then a vector containing
            -1 is returned
         */
        public static Vector<double> LCSS(Matrix<double> trajectory1, Matrix<double> trajectory2, int pointSpacing = -1, int pointDistance = 20, double errorMargin = 2){
            if(!TrajCheck(trajectory1, trajectory2)){
                return Vector<double>.Build.Dense(1, -1);
            }

            int dimensions = trajectory1.ColumnCount;
            int length1 = trajectory1.RowCount;
            int length2 = trajectory2.RowCount;

            if(length1 == 0 || length2 == 0){
                Console.WriteLine("At least one trajectory contains 0 points");
                return Vector<double>.Build.Dense(1, 0);
            }

            if(dimensions == 0){
                Console.WriteLine("Dimension is 0");
                return Vector<double>.Build.Dense(1, Math.Min(length1, length2));;
            }

            if(pointSpacing < 0){
                pointSpacing = Math.Max(length1, length2);
            }

            List<Vector<double>> translations = new List<Vector<double>>();
            //Calculating the subsets of translations
            for(int i = 0; i < dimensions; i++){
                translations.Add(TranslationSubset(trajectory1.Column(i), trajectory2.Column(i), pointSpacing, pointDistance));
            }

            //Storing the most optimal translations and similarity found so far
            int similarity = LCSSCalc(trajectory1, trajectory2, pointSpacing, pointDistance);
            Vector<double> optimalTrans = Vector<double>.Build.Dense(dimensions + 1, 0);
            List<double> similarityList = new List<double>(similarity);
            similarityList.AddRange(optimalTrans.ToArray());
            Vector<double> similarityVector = Vector<double>.Build.Dense(similarityList.ToArray());
            //Calculating how many translation possibilities are skipped for every one that is checked
            //using the error margin given
            double spacing = (double)translations[0].Count / (4.0 * (double)pointSpacing / errorMargin);

            if(spacing < 1){
                spacing = 1;
            }else if(spacing > (double)translations[0].Count / 2.0){
                spacing = (double)translations[0].Count / 2.0;
            }
            //Running the LCSS algorithm on each of the translations to be checked
            similarityVector = SimLoop(trajectory1, trajectory2, pointSpacing, pointDistance, (int)spacing, similarityVector, translations, dimensions, dimensions);

            return similarityVector;
        }

        /*Loops over and tests the trajectories using different translations in each dimension
        Args: 
            trajectory1: An m x n matrix containing trajectory1. Here m is the number 
                of points and n is the dimension of the points.
            trajectory2: traj2: A k x n matrix containing trajectory2. Here k is the number 
                of points and n is the dimension of the points.
                The two trajectories do not need to have the same number of points
            pointSpacing: An integer value of the maximum index difference between trajectory1 and
                trajectory2 allowed in the calculation. A negative values ets the spacing to unlimited
            pointDistance: A number representing the maximum distance in each dimension is allowed
                for points to be considered equivalent
            spacing: The integer spacing between each translation that will be tested
            similarity: A vector containing the current best similarity in each dimension
            dimensions: The number of dimensions used in the calculation
            dimLeft: The dimensions which have not been looped over yet
            currentTrans: A vector containing the current translation being tested
        Returns:
            Returns the current best LCSS value and the translations that created this as a vector
         */
        private static Vector<double> SimLoop(Matrix<double> trajectory1, Matrix<double> trajectory2, int pointSpacing, int pointDistance, int spacing, Vector<double> similarity, List<Vector<double>> translations, int dimensions, int dimLeft, Vector<double> currentTrans = null){
            if(currentTrans == null){
                currentTrans = Vector<double>.Build.Dense(dimensions, 0);
            }

            int thisDim = dimensions - dimLeft;

            double prevTrans = -1;
            //Testing each translation in this dimension
            for(int i = spacing - 1; i < translations[thisDim].Count; i += spacing){
                //The newest translation
                currentTrans[thisDim] = translations[thisDim][i];
                //Skipping translations which have already been checked
                if(currentTrans[thisDim] != prevTrans){
                    if(dimLeft > 1){
                        similarity = SimLoop(trajectory1, trajectory2, pointSpacing, pointDistance, spacing, similarity, translations, dimensions, dimLeft - 1, currentTrans);
                    }else{
                        //Running the LCSS algorithm on each of the translations to be checked
                        int newValue = LCSSCalc(trajectory1, trajectory2, pointSpacing, pointDistance, currentTrans);

                        //Keeping the similarity and translations if they're better than the previous best
                        if(newValue > similarity[0]){
                            similarity[0] = newValue;
                            for(int d = 0; d < dimensions; d++){
                                similarity[d + 1] = currentTrans[d];
                            }
                        }
                    }
                    prevTrans = currentTrans[thisDim];
                }
            }

            return similarity;
        }

        /*Calculates the LCSS of two trajectories using a set translation
        Args: 
            trajectory1: An m x n matrix containing trajectory1. Here m is the number 
                of points and n is the dimension of the points.
            trajectory2: traj2: A k x n matrix containing trajectory2. Here k is the number 
                of points and n is the dimension of the points.
                The two trajectories do not need to have the same number of points
            pointSpacing: An integer value of the maximum index difference between trajectory1 and
                trajectory2 allowed in the calculation. A negative values ets the spacing to unlimited
            pointDistance: A number representing the maximum distance in each dimension is allowed
                for points to be considered equivalent
            translations: A vector containing translations in each dimension to be applied to trajectory2
                in this calculation
        Returns:
            An integer representing the maximum LCSS value obtained using the variables provided. If a problem
            occurs then -1 is returned
         */
        public static int LCSSCalc(Matrix<double> trajectory1, Matrix<double> trajectory2, int pointSpacing = -1, int pointDistance = 20, Vector<double> translations = null){
            if(!TrajCheck(trajectory1, trajectory2)){
                return -1;
            }

            int dimensions = trajectory1.ColumnCount;
            int length1 = trajectory1.RowCount;
            int length2 = trajectory2.RowCount;

            if(translations == null){
                translations = Vector<double>.Build.Dense(dimensions, 0.0);
            }

            if(length1 == 0 || length2 == 0){
                Console.WriteLine("At least one trajectory contains 0 points");
                return 0;
            }

            if(dimensions == 0){
                Console.WriteLine("Dimension is 0");
                return Math.Min(length1, length2);
            }

            if(pointSpacing < 0){
                pointSpacing = Math.Max(length1, length2);
            }

            Matrix<double> distMatrix = Matrix<double>.Build.Dense(length1, length2, 0);
            int similarity = 0;

            for(int row = 0; row < length1; row++){
                //Calculating the relevant columns for each row
                int minCol = 0;
                int maxCol = length2 - 1;
                if(row > pointSpacing){
                    minCol = row - pointSpacing;
                }
                if(row < length2 - pointSpacing - 1){
                    maxCol = row + pointSpacing;
                }
                if(minCol <= maxCol){
                    for(int col = minCol; col <= maxCol; col++){
                        double newValue = 0;
                        double finalValue = 0;
                        //Calculating the new LCSS value for the current two points
                        //Checking the diagonal
                        if(row != 0 && col != 0){
                            finalValue = newValue = distMatrix[row - 1, col - 1];
                        }
                        //Checking below
                        if(row != 0){
                            double below = distMatrix[row - 1, col];
                            finalValue = Math.Max(below, finalValue);
                        }
                        //Checking the left
                        if(col != 0){
                            double before = distMatrix[row, col - 1];
                            finalValue = Math.Max(before, finalValue);
                        }
                        //Checking if the current points can increment the LCSS
                        if(finalValue < newValue + 1){
                            bool checkPoint = DistanceCheck(trajectory1.Row(row), trajectory2.Row(col) + translations, pointDistance);

                            if(checkPoint){
                                finalValue = ++newValue;
                            }
                        }
                        //Updating the distance matrix
                        distMatrix[row, col] = finalValue;
                        //Updating the similarity if a new maximum has been found
                        if(finalValue > similarity){
                            similarity = (int)finalValue;
                        }
                    }
                }          
            }
            return similarity;
        }

        /*Calculates the ratio of the longest common subsequence to the shortest trajectory
        Args: 
            trajectory1: An m x n matrix containing trajectory1. Here m is the number 
                of points and n is the dimension of the points.
            trajectory2: traj2: A k x n matrix containing trajectory2. Here k is the number 
                of points and n is the dimension of the points.
                The two trajectories do not need to have the same number of points
            pointSpacing: An integer value of the maximum index difference between trajectory1 and
                trajectory2 allowed in the calculation. A negative values ets the spacing to unlimited
            pointDistance: A number representing the maximum distance in each dimension is allowed
                for points to be considered equivalent
            errorMargin: A number used to scale the accuracy and speed of the calculation
        Returns:
            A number representing the maximum LCSS ratio obtained using the variables provided. If a problem occurs
            then -1 is returned
         */
        public static double LCSSRatio(Matrix<double> trajectory1, Matrix<double> trajectory2, int pointSpacing = -1, int pointDistance = 20, double errorMargin = 2){
            if(!TrajCheck(trajectory1, trajectory2)){
                return -1;
            }

            int dimensions = trajectory1.ColumnCount;
            int length1 = trajectory1.RowCount;
            int length2 = trajectory2.RowCount;

            if(length1 == 0 || length2 == 0){
                Console.WriteLine("At least one trajectory contains 0 points");
                return 0;
            }

            double length = Math.Min(length1, length2);

            return LCSS(trajectory1, trajectory2, pointSpacing, pointDistance, errorMargin)[0] / length;
        }

        /*Calculates the LCSS ratio between two trajectories using a set translation
        Args: 
            trajectory1: An m x n matrix containing trajectory1. Here m is the number 
                of points and n is the dimension of the points.
            trajectory2: traj2: A k x n matrix containing trajectory2. Here k is the number 
                of points and n is the dimension of the points.
                The two trajectories do not need to have the same number of points
            pointSpacing: An integer value of the maximum index difference between trajectory1 and
                trajectory2 allowed in the calculation. A negative values ets the spacing to unlimited
            pointDistance: A number representing the maximum distance in each dimension is allowed
                for points to be considered equivalent
            translations: A vector containing the translations in each dimension to be applied to
                trajectory2 in this calculation
        Returns:
            A number representing the maximum LCSS ratio obtained using the variables provided. If a 
            problem occurs then -1 is returned.
         */
        public static double LCSSRatioCalc(Matrix<double> trajectory1, Matrix<double> trajectory2, int pointSpacing = -1, int pointDistance = 20, Vector<double> translations = null){
            if(!TrajCheck(trajectory1, trajectory2)){
                return -1;
            }

            int dimensions = trajectory1.ColumnCount;
            int length1 = trajectory1.RowCount;
            int length2 = trajectory2.RowCount;

            if(translations == null){
                translations = Vector<double>.Build.Dense(dimensions, 0.0);
            }

            if(length1 == 0 || length2 == 0){
                Console.WriteLine("At least one trajectory contains 0 points");
                return 0;
            }

            double length = Math.Min(length1, length2);

            return (double)LCSSCalc(trajectory1, trajectory2, pointSpacing, pointDistance, translations) / length;
        }

        /*Calculates the frechet distance when one trajectory contains only a single point
        Args: 
            trajectory1: An m x n matrix containing trajectory1. Here m is the number 
                of points and n is the dimension of the points.
            trajectory2: traj2: A k x n matrix containing trajectory2. Here k is the number 
                of points and n is the dimension of the points.
                The two trajectories do not need to have the same number of points
        Returns:
            A number representing the Frechet distance is returned
         */
        private static double SinglePointCalc(Matrix<double> trajectory1, Matrix<double> trajectory2){
            int dimensions = trajectory1.ColumnCount;
            int length1 = trajectory1.RowCount;
            int length2 = trajectory2.RowCount;
            double leashSq = 1.0;

            if(length1 == 1){
                for(int point2 = 0; point2 < length2; point2++){
                    double newLeashSq = DistanceSq(trajectory1.Row(0), trajectory2.Row(point2));
                    leashSq = Math.Max(leashSq, newLeashSq);
                }
            }else if(length2 == 1){
                for(int point1 = 0; point1 < length1; point1++){
                    double newLeashSq = DistanceSq(trajectory1.Row(point1), trajectory2.Row(0));
                    leashSq = Math.Max(leashSq, newLeashSq);
                }
            }

            if(leashSq >= 0){
                return Math.Sqrt(leashSq);
            }else{
                Console.WriteLine("Error in single point trajectory calculation");
                return -1;
            }
        }

        /*Calculates the subsets of translations to be tested using the LCSS methods
        Args: 
            trajectory1: An m x n matrix containing trajectory1. Here m is the number 
                of points and n is the dimension of the points.
            trajectory2: traj2: A k x n matrix containing trajectory2. Here k is the number 
                of points and n is the dimension of the points.
                The two trajectories do not need to have the same number of points
            pointSpacing: An integer value of the maximum index difference between trajectory1 and
                trajectory2 allowed in the calculation. A negative values ets the spacing to unlimited
            pointDistance: A number representing the maximum distance in each dimension is allowed
                for points to be considered equivalent
        Retuns:
            A vector of numbers containing the translations calculated. This vector is sorted
            in ascending order
         */
        private static Vector<double> TranslationSubset(Vector<double> trajectory1, Vector<double> trajectory2, int pointSpacing, int pointDistance){
            int length1 = trajectory1.Count;
            int length2 = trajectory2.Count;

            //To turn into array, then into vector later
            List<double> translations = new List<double>();

            for(int row = 0; row < length1; row++){
                int minCol = 0;
                int maxCol = length2 - 1;

                if(row > pointSpacing){
                    minCol = row - pointSpacing;
                }

                if(row < length2 - pointSpacing - 1){
                    maxCol = row + pointSpacing;
                }

                if(minCol <= maxCol){
                    for(int col = minCol; col <= maxCol; col++){
                        translations.Add(trajectory1[row] - trajectory2[col] + pointDistance);
                        translations.Add(trajectory1[row] - trajectory2[col] - pointDistance);
                    }
                }
            }
            translations.Sort();

            return Vector<double>.Build.Dense(translations.ToArray());
        }

        /*Calculates a variation of trajectory2 by aligning its start and end points with those of trajectory1
        Args: 
            trajectory1: An m x n matrix containing trajectory1. Here m is the number 
                of points and n is the dimension of the points.
            trajectory2: traj2: A k x n matrix containing trajectory2. Here k is the number 
                of points and n is the dimension of the points.
                The two trajectories do not need to have the same number of points
        Returns:
            An m x n matrix containing the new variation of trajectory2 is returned.
            m is the number of points and n is the dimension of the points
         */
        public static Matrix<double> StartEndTranslate(Matrix<double> trajectory1, Matrix<double> trajectory2){
            if(!TrajCheck(trajectory1, trajectory2)){
                return null;
            }

            int dimensions = trajectory1.ColumnCount;
            int length1 = trajectory1.RowCount;
            int length2 = trajectory2.RowCount;

            if(length1 == 0 || length2 == 0){
                Console.WriteLine("A trajectory has no points");
                return trajectory2;
            }

            if(dimensions == 0){
                Console.WriteLine("Dimension is 0");
                return trajectory2;
            }

            Matrix<double> newTraj = Matrix<double>.Build.DenseOfMatrix(trajectory2);

            for(int i = 0; i < dimensions; i++){
                double diff1 = trajectory1[length1 - 1, i] - trajectory1[0, i];
                double diff2 = trajectory2[length2 - 1, i] - trajectory2[0, i];

                if(diff2 == 0){
                    Console.WriteLine("Equivalent start/end points in 1 dimension");
                }else{
                    for(int point = 0; point < length2; point++){
                        double pointDiff = trajectory2[point, i] - trajectory2[0, i];
                        newTraj[point, i] = (pointDiff / diff2) * diff1 + trajectory1[0, i];
                    }
                }
            }

            return newTraj;
        }

        /*Calculates the Frechet distance between two trajectories
        Args: 
            trajectory1: An m x n matrix containing trajectory1. Here m is the number 
                of points and n is the dimension of the points.
            trajectory2: traj2: A k x n matrix containing trajectory2. Here k is the number 
                of points and n is the dimension of the points.
                The two trajectories do not need to have the same number of points
            testLeash: A numeric leash value, which if positive, checks whether the leash
                can be used between the two trajectories. If this value is negative, then it is not
                used and the standard calculation is performed.
        Returns:
            A number representing the Frechet distance is returned. If a test leash is given, it
            is returned if the leash was successful, otherwise -1 is returned. If a problem occurs
            then -1 is returned.
         */
        public static double Frechet(Matrix<double> trajectory1, Matrix<double> trajectory2, double testLeash = -1.0){
            if(!TrajCheck(trajectory1, trajectory2)){
                return -1;
            }

            int dimensions = trajectory1.ColumnCount;
            int length1 = trajectory1.RowCount;
            int length2 = trajectory2.RowCount;

            if(length1 == 0 || length2 == 0){
                Console.WriteLine("At least one length is 0");
                return 0;
            }

            if(dimensions == 0){
                Console.WriteLine("The dimension is 0");
                return 0;
            }

            if(length1 == 1 || length2 == 1){
                double leash = SinglePointCalc(trajectory1, trajectory2);
                if(testLeash >= 0){
                    return Math.Min(testLeash, leash);
                }else{
                    return leash;
                }
            }

            Vector<double> dist1 = Vector<double>.Build.Dense(length1 - 1, 0);
            Vector<double> dist2 = Vector<double>.Build.Dense(length2 - 1, 0);

            for(int point = 0; point < length1 - 1; point++){
                double dist = DistanceSq(trajectory1.Row(point + 1), trajectory1.Row(point));
                dist1[point] = Math.Sqrt(dist);
            }

            for(int point = 0; point < length2 - 1; point++){
                double dist = DistanceSq(trajectory2.Row(point + 1), trajectory2.Row(point));
                dist2[point] = Math.Sqrt(dist);
            }

            double minLeashSq = DistanceSq(trajectory1.Row(0), trajectory2.Row(0));
            double endDistSq  = DistanceSq(trajectory1.Row(length1 - 1), trajectory2.Row(length2 - 1));
            minLeashSq = Math.Min(minLeashSq, endDistSq);

            List<double> leashList = new List<double>();
            leashList.Add(Math.Sqrt(minLeashSq));

            Matrix<double> distSq12 = Matrix.Build.Dense(length1, length2, 0);

            for(int point1 = 0; point1 < length1; point1++){
                for(int point2 = 0; point2 < length2; point2++){
                    double dist = DistanceSq(trajectory1.Row(point1), trajectory2.Row(point2));
                    distSq12[point1, point2] = dist;
                    //Adding in these leash possibilities because they are critical points
                    if(dist > minLeashSq){
                        leashList.Add(Math.Sqrt(dist));
                    }
                }
            }

            //If a testLeash is given
            if(testLeash >= 0){
                if(FrechetCheck(trajectory1, trajectory2, testLeash, dist1, dist2, distSq12)){
                    return testLeash;
                }else{
                    return -1;
                }
            }

            //Adding critical point leash possibilities to leash list
            for(int point1 = 0; point1 < length1 - 1; point1++){
                //Creating a unit vector in the direction of the next point from point1
                Vector<double> unitV1 = Vector<double>.Build.Dense(dimensions, 0);
                if(dist1[point1] != 0){
                    unitV1 = (trajectory1.Row(point1 + 1) - trajectory1.Row(point1)) / dist1[point1];
                }

                for(int point2 = 0; point2 < length2; point2++){
                    //Creating a vector from point1 to point2
                    Vector<double> vect12 = trajectory2.Row(point2) - trajectory1.Row(point1);
                    //Dot product finds how far from point1 the closest point on the line is
                    double pointDistance = unitV1.DotProduct(vect12);
                    //Square for easy calculation
                    double pointDistanceSq = pointDistance * pointDistance;
                    //The square of the distance between the line segment and the point
                    double shortDistance = distSq12[point1, point2] - pointDistanceSq;
                    double leashSq = 0;

                    if(pointDistance < 0){
                    }else if(pointDistance > dist1[point1]){
                        leashSq = distSq12[point1 + 1, point2];
                    }else{
                        leashSq = shortDistance;
                    }

                    //Adding the leash possibility to the list
                    if(leashSq > minLeashSq){
                        leashList.Add(Math.Sqrt(leashSq));
                    }
                }
            }

            for(int point2 = 0; point2 < length2 - 1; point2++){
                //Creating a unit vector in the direction of the next point from point1
                Vector<double> unitV1 = Vector<double>.Build.Dense(dimensions, 0);
                if(dist2[point2] != 0){
                    unitV1 = (trajectory2.Row(point2 + 1) - trajectory2.Row(point2)) / dist2[point2];
                }

                for(int point1 = 0; point1 < length1; point1++){
                    //Creating a vector from point1 to point2
                    Vector<double> vect12 = trajectory1.Row(point1) - trajectory2.Row(point2);
                    //Dot product finds how far from point1 the closest point on the line is
                    double pointDistance = unitV1.DotProduct(vect12);
                    //Square for easy calculation
                    double pointDistanceSq = pointDistance * pointDistance;
                    //The square of the distance between the line segment and the point
                    double shortDistance = distSq12[point1, point2] - pointDistanceSq;
                    double leashSq = 0;

                    if(pointDistance < 0){
                    }else if(pointDistance > dist2[point2]){
                        leashSq = distSq12[point1, point2 + 1];
                    }else{
                        leashSq = shortDistance;
                    }

                    //Adding the leash possibility to the list
                    if(leashSq > minLeashSq){
                        leashList.Add(Math.Sqrt(leashSq));
                    }
                }
            }

            //Calculating the critical points where new passages may open
            if(length1 > 3){
                for(int point2 = 0; point2 < length2 - 1; point2++){
                    //Creating a unit vector in the direction of the next point from point2
                    Vector<double> unitV2 = Vector<double>.Build.Dense(dimensions, 0);
                    if(dist2[point2] != 0){
                        unitV2 = (trajectory2.Row(point2 + 1) - trajectory2.Row(point2)) / dist2[point2];
                    }

                    for(int point1 = 1; point1 < length1 - 2; point1++){
                        //Creating a vector from point2 to point 1
                        Vector<double> vect21 = trajectory1.Row(point1) - trajectory2.Row(point2);
                        //Dot product finds how far from point2 the closest point on the line is
                        double pointDistance = unitV2.DotProduct(vect21);
                        if(pointDistance > 0){
                            //Square for easy calculation
                            double pointDistanceSq = pointDistance * pointDistance;
                            //The square of the distance between the line segment and the point
                            double shortDistance = distSq12[point1, point2] - pointDistanceSq;
                            //The second point where the passage opens up
                            for(int newPoint = point1 + 1; newPoint < length1 - 1; newPoint++){
                                //Creating new vector from point2 to newPoint
                                Vector<double> vect2new = trajectory1.Row(newPoint) - trajectory2.Row(point2);
                                //Dot product finds how far from point2 the closest point on the line is
                                double newPointDistance = unitV2.DotProduct(vect2new);
                                if(newPointDistance > pointDistance){
                                    double newPointDistSq = newPointDistance * newPointDistance;
                                    double newShortDistance = distSq12[newPoint, point2] - newPointDistSq;
                                    //The distance between the two closest points on the line
                                    double pointDiff = pointDistance - newPointDistance;
                                    //Finding the point where the passage opens
                                    double equivPoint = (pointDiff * pointDiff + shortDistance - newShortDistance) / (pointDiff * 2.0);

                                    if(equivPoint > 0 && equivPoint < dist2[point2]){
                                        double leashSq = newShortDistance + (equivPoint * equivPoint);
                                        if(leashSq > minLeashSq){
                                            leashList.Add(Math.Sqrt(leashSq));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if(length2 > 3){
                for(int point1 = 0; point1 < length1 - 1; point1++){
                    //Creating a unit vector in the direction of the next point from point2
                    Vector<double> unitV1 = Vector<double>.Build.Dense(dimensions, 0);
                    if(dist1[point1] != 0){
                        unitV1 = (trajectory2.Row(point1 + 1) - trajectory2.Row(point1)) / dist2[point1];
                    }

                    for(int point2 = 1; point2 < length2 - 2; point2++){
                        //Creating a vector from point1 to point2
                        Vector<double> vect12 = trajectory2.Row(point2) - trajectory1.Row(point1);
                        //Dot product finds how far from point2 the closest point on the line is
                        double pointDistance = unitV1.DotProduct(vect12);
                        if(pointDistance > 0){
                            //Square for easy calculation
                            double pointDistanceSq = pointDistance * pointDistance;
                            //The square of the distance between the line segment and the point
                            double shortDistance = distSq12[point1, point2] - pointDistanceSq;
                            //The second point where the passage opens up
                            for(int newPoint = point2 + 1; newPoint < length2 - 1; newPoint++){
                                //Creating new vector from point2 to newPoint
                                Vector<double> vect1new = trajectory2.Row(newPoint) - trajectory1.Row(point1);
                                //Dot product finds how far from point2 the closest point on the line is
                                double newPointDistance = unitV1.DotProduct(vect1new);
                                if(newPointDistance > pointDistance){
                                    double newPointDistSq = newPointDistance * newPointDistance;
                                    double newShortDistance = distSq12[point1, newPoint] - newPointDistSq;
                                    //The distance between the two closest points on the line
                                    double pointDiff = pointDistance - newPointDistance;
                                    //Finding the point where the passage opens
                                    double equivPoint = (pointDiff * pointDiff + shortDistance - newShortDistance) / (pointDiff * 2.0);

                                    if(equivPoint > 0 && equivPoint < dist1[point1]){
                                        double leashSq = newShortDistance + (equivPoint * equivPoint);
                                        if(leashSq > minLeashSq){
                                            leashList.Add(Math.Sqrt(leashSq));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            //Sort leash and remove duplicates
            leashList.Sort();
            List<double> uniqueLeash = new List<double>();
            uniqueLeash.Add(leashList[0]);
            double lastLeash = uniqueLeash[0];
            foreach(double item in leashList){
                if(lastLeash != item){
                    lastLeash = item;
                    uniqueLeash.Add(item);
                }
            }
            
            //Set up binary search for the list
            int startSearch = 0;
            int endSearch = uniqueLeash.Count - 1;
            //Making sure the largest leash is large enough
            if(FrechetCheck(trajectory1, trajectory2, uniqueLeash[endSearch], dist1, dist2, distSq12)){
                //Execute binary search
                while(startSearch < endSearch){
                    int current = (endSearch - startSearch / 2) * startSearch;
                    if(FrechetCheck(trajectory1, trajectory2, uniqueLeash[current], dist1, dist2, distSq12)){
                        endSearch = current;
                    }else{
                        startSearch = current + 1;
                    }
                }
                //Return the shortest leash for the trajectories
                return uniqueLeash[endSearch];
            }else{
                Console.WriteLine("Unable to find frechet distance");
                return -1;
            }
        }

        /*Checks whether a leash will work between two trajectories using the Frechet method
        Args: 
            trajectory1: An m x n matrix containing trajectory1. Here m is the number 
                of points and n is the dimension of the points.
            trajectory2: traj2: A k x n matrix containing trajectory2. Here k is the number 
                of points and n is the dimension of the points.
                The two trajectories do not need to have the same number of points
            leash: A leash value to be checked
            dist1: A vector containing the distance between each successive two points in trajectory1
            dist2: A vector containing the distance between each successive two points in trajectory2
            distSq12: A matrix containing the distance between each pair of two points where 1 point
                lies in trajectory1 and the other in trajectory2
        Returns:
            A boolean value, true if the leash is successful, false otherwise.
         */
        private static bool FrechetCheck(Matrix<double> trajectory1, Matrix<double> trajectory2, double leash, Vector<double> dist1, Vector<double> dist2, Matrix<double> distSq12){
            double leashSq = leash * leash;
            int dimensions = trajectory1.ColumnCount;
            int length1 = trajectory1.RowCount;
            int length2 = trajectory2.RowCount;

            double[, ,] left = new double[length1, length2 - 1, 2];
            double[, ,] bottom = new double[length1 - 1, length2, 2];
            double[, ,] newLeft = new double[length1, length2 - 1, 2];
            double[, ,] newBottom = new double[length1 - 1, length2, 2];

            if(leashSq < distSq12[0, 0] | leashSq < distSq12[length1 - 1, length2 - 2]){
                return false;
            }

            //Calculating the freespace of the first trajectory wrt the second
            for(int point1 = 0; point1 < length1 - 1; point1++){
                Vector<double> unitV1 = Vector<double>.Build.Dense(dimensions, 0);
                if(dist1[point1] != 0){
                    unitV1 = (trajectory1.Row(point1 + 1) - trajectory1.Row(point1) / dist1[point1]);
                }

                for(int point2 = 0; point2 < length2; point2++){
                    //Create vector from point 1 to point 2
                    Vector<double> vect12 = trajectory2.Row(point2) - trajectory1.Row(point1);
                    //Dot product finds how far from point1 the closest point on the line is
                    double pointDistance = unitV1.DotProduct(vect12);
                    //Square for easy calculation
                    double pointDistanceSq = pointDistance * pointDistance;
                    //Square of the distance between the line segment and the point
                    double shortDistance = distSq12[point1, point2] - pointDistanceSq;
                    //If some part of the current line can be used by the leash
                    if(shortDistance <= leashSq){
                        //Calculating the envelope alone the line
                        double envSize = Math.Sqrt(leashSq - shortDistance);
                        double envLow = pointDistance - envSize;
                        double envHigh = pointDistance + envSize;

                        //If whole line is in envelope
                        if(envHigh >= dist1[point1] && envLow <= 0){
                            bottom[point1, point2, 0] = 0;
                            bottom[point1, point2, 1] = 1;
                        }else if(envHigh >= 0 && envLow <= 0){
                            //If the start of the line is within the envelope
                            bottom[point1, point2, 0] = 0;
                            bottom[point1, point2, 1] = envHigh / dist1[point1];
                        }else if(envHigh >= dist1[point1] && envLow <= dist1[point1]){
                            //If the end of the line is within the envelope
                            bottom[point1, point2, 0] = envLow / dist1[point1];
                            bottom[point1, point2, 1] = 1;
                        }else if(envHigh >= 0 && envLow <= dist1[point1]){
                            //If the envelope is completely within the line
                            bottom[point1, point2, 0] = envLow / dist1[point1];
                            bottom[point1, point2, 1] = envHigh / dist1[point1];
                        }
                    }
                }
            }

            //Calculating the freespace of the second trajectory wrt the first
            for(int point2 = 0; point2 < length2 - 1; point2++){
                Vector<double> unitV1 = Vector<double>.Build.Dense(dimensions, 0);
                if(dist1[point2] != 0){
                    unitV1 = (trajectory2.Row(point2 + 1) - trajectory2.Row(point2) / dist2[point2]);
                }

                for(int point1 = 0; point1 < length1; point1++){
                    //Create vector from point 1 to point 2
                    Vector<double> vect12 = trajectory1.Row(point1) - trajectory2.Row(point2);
                    //Dot product finds how far from point1 the closest point on the line is
                    double pointDistance = unitV1.DotProduct(vect12);
                    //Square for easy calculation
                    double pointDistanceSq = pointDistance * pointDistance;
                    //Square of the distance between the line segment and the point
                    double shortDistance = distSq12[point1, point2] - pointDistanceSq;
                    //If some part of the current line can be used by the leash
                    if(shortDistance <= leashSq){
                        //Calculating the envelope alone the line
                        double envSize = Math.Sqrt(leashSq - shortDistance);
                        double envLow = pointDistance - envSize;
                        double envHigh = pointDistance + envSize;

                        //If whole line is in envelope
                        if(envHigh >= dist2[point2] && envLow <= 0){
                            left[point1, point2, 0] = 0;
                            left[point1, point2, 1] = 1;
                        }else if(envHigh >= 0 && envLow <= 0){
                            //If the start of the line is within the envelope
                            left[point1, point2, 0] = 0;
                            left[point1, point2, 1] = envHigh / dist2[point2];
                        }else if(envHigh >= dist2[point2] && envLow <= dist2[point2]){
                            //If the end of the line is within the envelope
                            left[point1, point2, 0] = envLow / dist2[point2];
                            left[point1, point2, 1] = 1;
                        }else if(envHigh >= 0 && envLow <= dist2[point2]){
                            //If the envelope is completely within the line
                            left[point1, point2, 0] = envLow / dist2[point2];
                            left[point1, point2, 1] = envHigh / dist2[point2];
                        }
                    }
                }
            }
            //Set up new arrays to find the monotone freespace
            newLeft[0, 0, 0] = left[0, 0, 0];
            newLeft[0, 0, 1] = left[0, 0, 1];
            newBottom[0, 0, 0] = bottom[0, 0, 0];
            newBottom[0, 0, 1] = bottom[0, 0, 1];

            //Setting the first line of the new left array
            if(length2 > 2){
                for(int point2 = 1; point2 < length2 - 1; point2++){
                    if(newLeft[0, point2 - 1, 1] == 1){
                        newLeft[0, point2, 0] = left[0, point2, 0];
                        newLeft[0, point2, 1] = left[0, point2, 1];
                    }
                }
            }

            //Setting the first line of the new bottom array
            if(length1 > 2){
                for(int point1 = 1; point1 < length1 - 1; point1++){
                    if(newBottom[point1 - 1, 0, 1] == 1){
                        newBottom[point1, 0, 0] = bottom[point1, 0, 0];
                        newBottom[point1, 0, 1] = bottom[point1, 0, 1];
                    }
                }
            }

            //Calculating the monotone freespace
            for(int point1 = 0; point1 < length1; point1++){
                for(int point2 = 0; point2 < length2; point2++){
                    if(point1 != length1 - 1 && point2 != 0){
                        //If the area is allowable from the freespace
                        if(bottom[point1, point2, 0] > -0.1){
                            //If the new area can be entered from the left.
                            if(newLeft[point1, point2 - 1, 0] > -0.1){
                                //Setting up the montone freespace for these points.
                                newBottom[point1, point2, 0] = bottom[point1, point2, 0];
                                newBottom[point1, point2, 1] = bottom[point1, point2, 1];
                            }else if(newBottom[point1, point2 - 1, 0] > -0.1){
                                if(newBottom[point1, point2 - 1, 0] <= bottom[point1, point2, 0]){
                                    newBottom[point1, point2, 0] = bottom[point1, point2, 0];
                                    newBottom[point1, point2, 1] = bottom[point1, point2, 1];
                                }else if(newBottom[point1, point2 - 1, 0] <= bottom[point1, point2, 1]){
                                    newBottom[point1, point2, 0] = bottom[point1, point2 - 1, 0];
                                    newBottom[point1, point2, 1] = bottom[point1, point2, 1];
                                }
                            }
                        }
                    }
                    if(point2 != length2 - 1 && point1 != 0){
                        //If the area is allowable from the freespace
                        if(left[point1, point2, 0] > -0.1){
                            //If the new area can be entered from the bottom.
                            if(newBottom[point1 - 1, point2, 0] > -0.1){
                                //Setting up the montone freespace for these points.
                                newLeft[point1, point2, 0] = left[point1, point2, 0];
                                newLeft[point1, point2, 1] = left[point1, point2, 1];
                            }else if(newLeft[point1 - 1, point2, 0] > -0.1){
                                if(newLeft[point1 - 1, point2, 0] <= left[point1, point2, 0]){
                                    newLeft[point1, point2, 0] = left[point1, point2, 0];
                                    newLeft[point1, point2, 1] = left[point1, point2, 1];
                                }else if(newLeft[point1 - 1, point2, 0] <= left[point1, point2, 1]){
                                    newLeft[point1, point2, 0] = left[point1 - 1, point2 - 1, 0];
                                    newLeft[point1, point2, 1] = left[point1, point2, 1];
                                }
                            }
                        }
                    }
                }
            }
            //If the monotone freespace reaches the final point then the leash is successful
            if(newLeft[length1 - 1, length2 - 2, 1] == 1 || newBottom[length1 - 2, length2 - 1, 1] == 1){
                return true;
            }else{
                return false;
            }
        }
    }    
}