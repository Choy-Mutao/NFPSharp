using NFPSharp;
using QuickSVG;

namespace NFPTest
{
    public class Tests
    {
        private void DeleteFile(string path)
        {
            try
            {
                if (File.Exists(path))
                {
                    File.Delete(path);
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Error when delete file : {ex.Message}");
            }
        }

        [SetUp]
        public void Setup()
        {
        }

        [Test]
        public void Test_QuickSvg()
        {
            double[][] points =
            [
                [100, 10],
                [0, 0],
                [100, 90],
                [100, 10],
            ];

            QuickSVGWriter writer = new QuickSVGWriter();
            QuickSVGUtil.AddSolution(writer, new PathD(points), true);
            string filename = @"..\..\..\Test_QuickSvg.svg";
            QuickSVGUtil.SaveToFile(writer, filename, FillRule.NonZero, 800, 600, 10);
            var p = QuickSVGUtil.OpenFileWithDefaultApp(filename);
            DeleteFile(filename);
            Assert.IsTrue(true);
        }

        /// <summary>
        /// !!! A Tough Method to detect if a point on the segment by (A, B)
        /// This code from SVGNest projet, to test the new method result
        /// </summary>
        /// <returns></returns>
        private bool OnSegment(NFPPoint A, NFPPoint B, NFPPoint p, double tolerance = 1e-9)
        {

            // vertical line
            if (MathUtil.AlmostEqual(A.x, B.x, tolerance) && MathUtil.AlmostEqual(p.x, A.x, tolerance))
            {
                if (!MathUtil.AlmostEqual(p.y, B.y, tolerance) && !MathUtil.AlmostEqual(p.y, A.y, tolerance) && p.y < Math.Max(B.y, A.y) && p.y > Math.Min(B.y, A.y))
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }

            // horizontal line
            if (MathUtil.AlmostEqual(A.y, B.y, tolerance) && MathUtil.AlmostEqual(p.y, A.y, tolerance))
            {
                if (!MathUtil.AlmostEqual(p.x, B.x, tolerance) && !MathUtil.AlmostEqual(p.x, A.x, tolerance) && p.x < Math.Max(B.x, A.x) && p.x > Math.Min(B.x, A.x))
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }

            //range check, !!! too tight !!!
            if ((p.x < A.x && p.x < B.x) || (p.x > A.x && p.x > B.x) || (p.y < A.y && p.y < B.y) || (p.y > A.y && p.y > B.y))
            {
                return false;
            }


            // exclude end points
            if ((MathUtil.AlmostEqual(p.x, A.x, tolerance) && MathUtil.AlmostEqual(p.y, A.y, tolerance)) || (MathUtil.AlmostEqual(p.x, B.x, tolerance) && MathUtil.AlmostEqual(p.y, B.y, tolerance)))
            {
                return false;
            }

            var cross = (p.y - A.y) * (B.x - A.x) - (p.x - A.x) * (B.y - A.y);

            if (Math.Abs(cross) > tolerance)
            {
                return false;
            }

            var dot = (p.x - A.x) * (B.x - A.x) + (p.y - A.y) * (B.y - A.y);

            if (dot < 0 || MathUtil.AlmostEqual(dot, 0, tolerance))
            {
                return false;
            }

            var len2 = (B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y);

            if (dot > len2 || MathUtil.AlmostEqual(dot, len2, tolerance))
            {
                return false;
            }

            return true;
        }

        [Test]
        public void Test_OnSegment()
        {
            double TOL = 1e-9;

            NFPPoint A = new NFPPoint(-5, -2), B = new NFPPoint(2, -2);
            NFPPoint p = new NFPPoint(-1, -2);

            Assert.IsTrue(OnSegment(A, B, p, TOL) == true);
            Assert.IsTrue(MathUtil.OnSegment(A, B, p, TOL) == true);

            p = new NFPPoint(-5, -2);
            Assert.IsTrue(OnSegment(A, B, p, TOL) == false);
            Assert.IsTrue(MathUtil.OnSegment(A, B, p, TOL) == false);

            p = new NFPPoint(2, -2);
            Assert.IsTrue(OnSegment(A, B, p, TOL) == false);
            Assert.IsTrue(MathUtil.OnSegment(A, B, p, TOL) == false);

            p = new NFPPoint(-1.35427, -2.00054);
            Assert.IsTrue(OnSegment(A, B, p, TOL) == false);
            Assert.IsTrue(MathUtil.OnSegment(A, B, p, TOL) == false);

            Assert.IsTrue(OnSegment(A, B, p, 1) == true);
            Assert.IsTrue(MathUtil.OnSegment(A, B, p, 1e-7) == true);

            Assert.IsTrue(OnSegment(A, B, p, 1e-8) == false);
            Assert.IsTrue(MathUtil.OnSegment(A, B, p, 1e-8) == false);
        }

        [Test]
        public void Test_PointDistance()
        {
            NFPPoint p, s1, s2;

            // horizontal line
            p = new NFPPoint(0, 0);
            s1 = new NFPPoint(-10, 3);
            s2 = new NFPPoint(10, 3);

            NFPPoint point;

            Assert.IsTrue(NFPUtils.PointDistance(p, s1, s2, new NFPVector(0, 1)) == 3);
            Assert.IsTrue(NFPUtils.PointDistance(p, s1, s2, new NFPVector(1, 1)) == 3 * Math.Sqrt(2));
            Assert.IsTrue(NFPUtils.PointDistance(p, s1, s2, new NFPVector(-1, 1)) == 3 * Math.Sqrt(2));

            Assert.IsTrue(NFPUtils.PointDistance(p, s1, s2, new NFPVector(0, -1)) == -3);
            Assert.IsTrue(NFPUtils.PointDistance(p, s1, s2, new NFPVector(-1, -1)) == -3 * Math.Sqrt(2));
            Assert.IsTrue(NFPUtils.PointDistance(p, s1, s2, new NFPVector(1, -1)) == -3 * Math.Sqrt(2));

            Assert.IsTrue(NFPUtils.PointDistance(p, s1, s2, new NFPVector(-10, 3)) == null);
            Assert.IsTrue(NFPUtils.PointDistance(p, s1, s2, new NFPVector(10, 3)) == null);
        }

        [Test]
        public void Test_SegmentDistance()
        {
            NFPPoint A = new NFPPoint(10, 0);
            NFPPoint B = new NFPPoint(10, 10);

            NFPPoint E = new NFPPoint(15, 0);
            NFPPoint F = new NFPPoint(15, 10);

            var d = NFPUtils.SegmentDistance(A, B, E, F, new NFPVector(1, 0));
            Assert.IsTrue(d != null);
        }

        [Test]
        public void Test_PolygonSlideDistance()
        {
            double[][] a = [[0, 0], [10, 0], [10, 10], [0, 10]];
            double[][] b = [[15, 0], [20, 0], [20, 10], [15, 10]];

            NFPPolygon A = new NFPPolygon(a);
            NFPPolygon B = new NFPPolygon(b);

            var d = NFPUtils.PolygonSlideDistance(A, B, new NFPVector(1, 0), false);

            Assert.IsTrue(d != null);
            Assert.IsTrue(d!.Value == 5);

            QuickSVGWriter writer = new QuickSVGWriter();
            QuickSVGUtil.AddSolution(writer, new PathD(a), true);
            QuickSVGUtil.AddSolution(writer, new PathD(b), true);
            string filename = @"..\..\..\Test_PolygonSlideDistance.svg";
            QuickSVGUtil.SaveToFile(writer, filename, FillRule.NonZero, 800, 600, 10);
            var p = QuickSVGUtil.OpenFileWithDefaultApp(filename);
            DeleteFile(filename);
            Assert.IsTrue(true);
        }

        [Test]
        public void Test_PointInPolygon()
        {
            double[][] a = [[0, 0], [10, 0], [10, 10], [0, 10]];
            double[][] b = [[15, 0], [20, 0], [20, 5], [15, 5]];

            NFPPolygon A = new NFPPolygon(a), B = new NFPPolygon(b);
            //A.Add(A.First()); B.Add(B.First());

            // drawing test
            QuickSVGWriter writer = new QuickSVGWriter();
            QuickSVGUtil.AddSubject(writer, new PathD(a));
            QuickSVGUtil.AddSubject(writer, new PathD(b));

            int in_count = 0;
            for (int i = 0; i < A.Count; i++)
            {
                if (!A[i].marked)
                {
                    A[i].marked = true;
                    for (int j = 0; j < B.Count; j++)
                    {
                        // 1. step 
                        B.offsetx = A[i].x - B[j].x;
                        B.offsety = A[i].y - B[j].y;

                        double[][] after_path = new double[B.Count][];
                        bool? Binside = null;
                        for(int k = 0; k < B.Count; k++)
                        {
                            var inpoly = NFPUtils.PointInPolygon(new NFPPoint(B[k].x + B.offsetx, B[k].y + B.offsety), A, 1e-9);
                            after_path[k] = [B[k].x + B.offsetx, B[k].y + B.offsety];

                            if (inpoly is not null && inpoly.Value == true) in_count++;
                        }
                        QuickSVGUtil.AddSolution(writer, new PathD(after_path), true);
                    }
                }
            }

            string filename = @"..\..\..\Test_PointInPolygon.svg";
            QuickSVGUtil.SaveToFile(writer, filename, FillRule.NonZero, 800, 600, 10);
            var p = QuickSVGUtil.OpenFileWithDefaultApp(filename);

            Assert.That(in_count, Is.EqualTo(4));
        }

        [Test]
        public void Test_PolygonProjectionDistance()
        {
            double[][] a = [[0, 0], [10, 0], [10, 10], [0, 10]];
            double[][] b = [[15, 0], [20, 0], [20, 5], [15, 5]];

            NFPPolygon A = new NFPPolygon(a);
            NFPPolygon B = new NFPPolygon(b);

            NFPVector dir = new NFPVector(1, 0.2);

            var result = NFPUtils.PolygonProjectionDistance(A, B, dir, 1e-9);
            if (result is not null)
            {
                NFPVector vec = dir;
                vec.SetLength(result.Value);

                double[][] after_path = B.Select(pnt => new double[] { pnt.x + vec.x, pnt.y + vec.y }).ToArray();

                // drawing test
                QuickSVGWriter writer = new QuickSVGWriter();
                QuickSVGUtil.AddSubject(writer, new PathD(a));
                QuickSVGUtil.AddSubject(writer, new PathD(b));
                QuickSVGUtil.AddSolution(writer, new PathD(after_path), true);
                string filename = @"..\..\..\Test_NoFitPolygon.svg";
                QuickSVGUtil.SaveToFile(writer, filename, FillRule.NonZero, 800, 600, 10);
                var p = QuickSVGUtil.OpenFileWithDefaultApp(filename);
            }

            Assert.That(result, Is.Not.Null);
        }


        [Test]
        public void Test_SearchStartPoint()
        {
            double[][] a = [[0, 0], [10, 0], [10, 10], [0, 10]];
            double[][] b = [[15, 0], [20, 0], [20, 5], [15, 5]];

            NFPPolygon A = new NFPPolygon(a);
            NFPPolygon B = new NFPPolygon(b);
            var result = NFPUtils.SearchStartPoint(A, B, true);

            if (result is not null)
            {
                var vec = result;
                double[][] after_path = B.Select(pnt => new double[] { pnt.x + vec.x, pnt.y + vec.y }).ToArray();

                // drawing test
                QuickSVGWriter writer = new QuickSVGWriter();
                QuickSVGUtil.AddSubject(writer, new PathD(a));
                QuickSVGUtil.AddSubject(writer, new PathD(b));
                QuickSVGUtil.AddSolution(writer, new PathD(after_path), true);
                string filename = @"..\..\..\Test_SearchStartPoint.svg";
                QuickSVGUtil.SaveToFile(writer, filename, FillRule.NonZero, 800, 600, 10);
                var p = QuickSVGUtil.OpenFileWithDefaultApp(filename);
            }

            Assert.Multiple(() =>
            {
                Assert.That(A, Has.Count.EqualTo(a.Length));
                Assert.That(result is not null, Is.True);
            });
        }

        [Test]
        public void Test_NoFitPolygon()
        {
            double[][] a = [[0, 0], [10, 0], [10, 10], [0, 10]];
            double[][] b = [[15, 0], [20, 0], [20, 5], [15, 5]];

            NFPPolygon A = new NFPPolygon(a);
            NFPPolygon B = new NFPPolygon(b);

            var result = NFPUtils.NoFitPolygon(A, B, true, true);
            var paths = result?.Select(item => new PathD(item.Select(pnt => new double[] { pnt.x, pnt.y })));

            QuickSVGWriter writer = new QuickSVGWriter();
            QuickSVGUtil.AddSubject(writer, new PathD(a));
            QuickSVGUtil.AddSubject(writer, new PathD(b));
            QuickSVGUtil.AddSolution(writer, new PathsD(paths!), true);
            string filename = @"..\..\..\Test_NoFitPolygon.svg";
            QuickSVGUtil.SaveToFile(writer, filename, FillRule.NonZero, 800, 600, 10);
            var p = QuickSVGUtil.OpenFileWithDefaultApp(filename);
            Assert.That(true, Is.True);
            //DeleteFile(filename);
        }
    }
}