using System;

namespace g4
{
    public struct Triangle3d
    {
        public Vector3d V0, V1, V2;

        public Triangle3d(Vector3d v0, Vector3d v1, Vector3d v2)
        {
            V0 = v0; V1 = v1; V2 = v2;
        }

        public Vector3d this[int key]
        {
            get { return (key == 0) ? V0 : (key == 1) ? V1 : V2; }
            set { if (key == 0) V0 = value; else if (key == 1) V1 = value; else V2 = value; }
        }

        public Vector3d Normal {
            get { return MathUtil.Normal(ref V0, ref V1, ref V2); }
        }
        public double Area {
            get { return MathUtil.Area(ref V0, ref V1, ref V2); }
        }
        public double AspectRatio {
            get { return MathUtil.AspectRatio(ref V0, ref V1, ref V2); }
        }

        public Vector3d PointAt(double bary0, double bary1, double bary2)
        {
            return bary0 * V0 + bary1 * V1 + bary2 * V2;
        }
        public Vector3d PointAt(Vector3d bary)
        {
            return bary.x* V0 + bary.y* V1 + bary.z* V2;
        }

        public Vector3d BarycentricCoords(Vector3d point)
        {
            return MathUtil.BarycentricCoords(point, V0, V1, V2);
        }

        public Vector3d GetCentroid()
        {
            return new Vector3d(
                (V0.x + V1.x + V2.x) / 3.0,
                (V0.y + V1.y + V2.y) / 3.0,
                (V0.z + V1.z + V2.z) / 3.0
            );
        }

        public bool IsPointInTriangle(Vector3d point, double epsilon = 1e-06)
        {
            return IsPointInTriangleBaryCoords(BarycentricCoords(point), epsilon);
        }

        private static bool IsPointInTriangleBaryCoords(Vector3d baryCoords, double epsilon = 1e-06)
        {
            // Check if all barycenter coordinates are between 0 and 1, inclusive
            var condition1 = baryCoords.x >= -epsilon && baryCoords.x <= 1.0 + epsilon;
            var condition2 = baryCoords.y >= -epsilon && baryCoords.y <= 1.0 + epsilon;
            var condition3 = baryCoords.z >= -epsilon && baryCoords.z <= 1.0 + epsilon;

            // Check if the sum of barycenter coordinates is approximately 1
            var sumIsOne = Math.Abs(baryCoords.x + baryCoords.y + baryCoords.z - 1.0) < epsilon;

            return condition1 && condition2 && condition3 && sumIsOne;
        }

        public bool IsDegenerate(float epsilon = 1e-4f)
        {
            // Check if any two vertices are coincident
            if (V0.Distance(V1) < epsilon ||
                V1.Distance(V2) < epsilon ||
                V0.Distance(V2) < epsilon)
            {
                return true;
            }

            // Check if the vertices are collinear
            var e1 = V1 - V0;
            var e2 = V2 - V0;

            // If the cross product is (almost) zero, the points are collinear
            return e1.Cross(e2).Length <= epsilon;
        }

        // conversion operators
        public static implicit operator Triangle3d(Triangle3f v)
        {
            return new Triangle3d(v.V0, v.V1, v.V2);
        }
        public static explicit operator Triangle3f(Triangle3d v)
        {
            return new Triangle3f((Vector3f)v.V0, (Vector3f)v.V1, (Vector3f)v.V2);
        }
    }



    public struct Triangle3f
    {
        public Vector3f V0, V1, V2;

        public Triangle3f(Vector3f v0, Vector3f v1, Vector3f v2)
        {
            V0 = v0; V1 = v1; V2 = v2;
        }

        public Vector3f this[int key]
        {
            get { return (key == 0) ? V0 : (key == 1) ? V1 : V2; }
            set { if (key == 0) V0 = value; else if (key == 1) V1 = value; else V2 = value; }
        }


        public Vector3f PointAt(float bary0, float bary1, float bary2)
        {
            return bary0 * V0 + bary1 * V1 + bary2 * V2;
        }
        public Vector3f PointAt(Vector3f bary)
        {
            return bary.x * V0 + bary.y * V1 + bary.z * V2;
        }

        public Vector3f BarycentricCoords(Vector3f point)
        {
            return (Vector3f)MathUtil.BarycentricCoords(point, V0, V1, V2);
        }

        public bool IsDegenerate(float epsilon = 1e-4f)
        {
            // Check if any two vertices are coincident
            if (V0.Distance(V1) < epsilon ||
                V1.Distance(V2) < epsilon ||
                V0.Distance(V2) < epsilon)
            {
                return true;
            }

            // Check if the vertices are collinear
            var e1 = V1 - V0;
            var e2 = V2 - V0;

            // If the cross product is (almost) zero, the points are collinear
            return e1.Cross(e2).Length <= epsilon;
        }
    }

}
