using NFPSharp;

namespace NFPTest
{
    public class Tests
    {
        [SetUp]
        public void Setup()
        {
        }

        [Test]
        public void Test1()
        {
            Assert.Pass();
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
    }
}