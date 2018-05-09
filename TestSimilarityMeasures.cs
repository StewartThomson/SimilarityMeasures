using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Collections.Generic;
using System.IO;
using CsvHelper;

namespace SimilarityMeasures
{
    public class coords{
        public double x {get; set;}
        public double y {get; set;}
    }
    public class TestSimilarityMeasures
    {
        public static void Main(){
            StreamReader reader = File.OpenText(@"E:\Unity\SimilarityMeasuresTestingUnity\Assets\CSV\curve1.csv");
            var csv = new CsvReader(reader);
            var records = csv.GetRecords<coords>();
            List<double[]> rows = new List<double[]>();

            foreach(var record in records){
                double[] row = new double[2];
                row[0] = record.x;
                row[1] = record.y;
                rows.Add(row);
            }
            Matrix<double> path1 = Matrix.Build.DenseOfRows(rows);
            Console.WriteLine(path1);
            reader.Close();

            reader = File.OpenText(@"E:\Unity\SimilarityMeasuresTestingUnity\Assets\CSV\curve2.csv");
            csv = new CsvReader(reader);
            records = csv.GetRecords<coords>();
            rows.Clear();
            
            foreach(var record in records){
                double[] row = new double[2];
                row[0] = record.x;
                row[1] = record.y;
                rows.Add(row);
            }
            Matrix<double> path2 = Matrix.Build.DenseOfRows(rows);
            reader.Close();
            Console.WriteLine(path2);
            
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
            Console.WriteLine("\nTESTING StartEndTranslate:\n" + SimilarityMeasures.StartEndTranslate(path1, path2));
            Console.WriteLine("\nTESTING Frechet:\n" + SimilarityMeasures.Frechet(path1, path2));
            Console.WriteLine("\nTESTING Frechcet with leash length 5:\n" + SimilarityMeasures.Frechet(path1, path2, 5));

            /*double[,] arr1 = {{0,1.6,17.85,3.25,27.3},{0,19.2,6.8,14.44,12.6}};
            path1 = Matrix.Build.DenseOfArray(arr1).Transpose();
            double[,] arr2 = {{0,1,2,3,4,5,11.3,40},{0,1,2,3,4,5,6.7,4}};
            path2 = Matrix.Build.DenseOfArray(arr2).Transpose();
            Console.WriteLine("\nTESTING Frechet:\n" + SimilarityMeasures.Frechet(path1, path2));
            Console.WriteLine("\nTESTING Frechcet with leash length 8:\n" + SimilarityMeasures.Frechet(path1, path2, 8));
            Console.WriteLine("\nTESTING StartEndTranslate:\n" + SimilarityMeasures.StartEndTranslate(path2, path1));
            Console.WriteLine("\nTESTING LCSS:\n" + SimilarityMeasures.LCSS(path1, path2, 2, 2, 0.5));*/
        }
    }
}