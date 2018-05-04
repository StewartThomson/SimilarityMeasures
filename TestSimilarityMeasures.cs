using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Collections.Generic;

namespace SimilarityMeasures
{
    public class TestSimilarityMeasures
    {
        public static void Main(){
            double[,] marr1 = {{0,1,2,3},{0,1,2,3}};
            Matrix<double> path1 = Matrix.Build.DenseOfArray(marr1).Transpose();
            double[,] marr2 = {{0,1,2,3},{4,5,6,7}};
            Matrix<double> path2 = Matrix.Build.DenseOfArray(marr2).Transpose();

            double[] varr1 = {0,2,4,6};
            Vector<double> point1 = Vector.Build.DenseOfArray(varr1);
            double[] varr2 = {0,1,2,3};
            Vector<double> point2 = Vector.Build.DenseOfArray(varr2);

            Console.WriteLine("\nTESTING AveTranslate:\n" + SimilarityMeasures.AveTranslate(path1, path2));
            Console.WriteLine("\nTESTING DistanceCheck with distance 1:\n" + SimilarityMeasures.DistanceCheck(point1, point2, 1));
            Console.WriteLine("\nTESTING DistanceCheck with distance 3:\n" + SimilarityMeasures.DistanceCheck(point1, point2, 3));
            Console.WriteLine("\nTESTING DistanceSq:\n" + SimilarityMeasures.DistanceSq(point1, point2));
            Console.WriteLine("\nTESTING DWT:\n" + SimilarityMeasures.DynamicTimeWarping(path1, path2));
            Console.WriteLine("\nTESTING EditDist:\n" + SimilarityMeasures.EditDist(path1, path2, 2));
            Console.WriteLine("\nTESTING LCSS:\n" + SimilarityMeasures.LCSS(path1, path2, 2, 2, 0.5));
        }
    }
}