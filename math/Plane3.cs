using System;

namespace g4
{
    //3D plane, based on WildMagic5 Wm5Plane3 class

    public struct Plane3d
    {
        public Vector3d Normal;
        public double Constant;

        public Plane3d(Vector3d normal, double constant)
        {
            Normal = normal;
            Constant = constant;
        }

         // N is specified, c = Dot(N,P) where P is a point on the plane.
        public Plane3d(Vector3d normal, Vector3d point)
        {
            Normal = normal;
            Constant = Normal.Dot(point);
        }

        // N = Cross(P1-P0,P2-P0)/Length(Cross(P1-P0,P2-P0)), c = Dot(N,P0) where
        // P0, P1, P2 are points on the plane.
        public Plane3d(Vector3d p0, Vector3d p1, Vector3d p2)
        {
            Vector3d edge1 = p1 - p0;
            Vector3d edge2 = p2 - p0;
            Normal = edge1.UnitCross(edge2);
            Constant = Normal.Dot(p0);
        }


        // Compute d = Dot(N,P)-c where N is the plane normal and c is the plane
        // constant.  This is a signed distance.  The sign of the return value is
        // positive if the point is on the positive side of the plane, negative if
        // the point is on the negative side, and zero if the point is on the
        // plane.
        public double DistanceTo(Vector3d p)
        {
            return Normal.Dot(p) - Constant;
        }

        // The "positive side" of the plane is the half space to which the plane
        // normal points.  The "negative side" is the other half space.  The
        // function returns +1 when P is on the positive side, -1 when P is on the
        // the negative side, or 0 when P is on the plane.
        public int WhichSide (Vector3d p)
        {
            double distance = DistanceTo(p);
            if (distance < 0)
                return -1;
            else if (distance > 0)
                return +1;
            else
                return 0;
        }

        /// <summary>
        /// Finds an arbitrary point on a plane defined by a normal and a constant (ax + by + cz - d = 0).
        /// </summary>
        /// <param name="epsilon">Precision.</param>
        /// <returns>A point that lies on the plane.</returns>
        public Vector3d GetAnyPoint(double epsilon = 1e-06)
        {
            // Extract components of the normal vector
            var a = Normal.x;
            var b = Normal.y;
            var c = Normal.z;

            // Choose arbitrary values for x and y
            var x = 0.0;
            var y = 0.0;
            var z = 0.0;

            // Solve for z using the plane equation
            if (Math.Abs(c) > epsilon)
            {
                z = -(a * x + b * y - Constant) / c;
            }
            else if (Math.Abs(b) > epsilon)
            {
                y = -(a * x - Constant) / b; // Solve for y if c is zero
            }
            else
            {
                x = Constant / a; // Solve for x if c and b are zero
            }

            return new Vector3d(x, y, z);
        }
    }




    public struct Plane3f
    {
        public Vector3f Normal;
        public float Constant;

        public Plane3f(Vector3f normal, float constant)
        {
            Normal = normal;
            Constant = constant;
        }

         // N is specified, c = Dot(N,P) where P is a point on the plane.
        public Plane3f(Vector3f normal, Vector3f point)
        {
            Normal = normal;
            Constant = Normal.Dot(point);
        }

        // N = Cross(P1-P0,P2-P0)/Length(Cross(P1-P0,P2-P0)), c = Dot(N,P0) where
        // P0, P1, P2 are points on the plane.
        public Plane3f(Vector3f p0, Vector3f p1, Vector3f p2)
        {
            Vector3f edge1 = p1 - p0;
            Vector3f edge2 = p2 - p0;
            Normal = edge1.UnitCross(edge2);
            Constant = Normal.Dot(p0);
        }


        // Compute d = Dot(N,P)-c where N is the plane normal and c is the plane
        // constant.  This is a signed distance.  The sign of the return value is
        // positive if the point is on the positive side of the plane, negative if
        // the point is on the negative side, and zero if the point is on the
        // plane.
        public float DistanceTo(Vector3f p)
        {
            return Normal.Dot(p) - Constant;
        }

        // The "positive side" of the plane is the half space to which the plane
        // normal points.  The "negative side" is the other half space.  The
        // function returns +1 when P is on the positive side, -1 when P is on the
        // the negative side, or 0 when P is on the plane.
        public int WhichSide (Vector3f p)
        {
            float distance = DistanceTo(p);
            if (distance < 0)
                return -1;
            else if (distance > 0)
                return +1;
            else
                return 0;
        }

    }


}
