using System;

namespace g4.intersection
{
    public static class IntrRay3Plane3
    {
        /// <summary>
        /// Finds the intersection point between a ray and a plane.
        /// </summary>
        /// <param name="ray">The ray defined by an origin and direction.</param>
        /// <param name="plane">The plane defined by a point and a normal.</param>
        /// <param name="epsilon">Precision.</param>
        /// <returns>
        /// A Point3D representing the intersection point, or null if there is no intersection.
        /// </returns>
        public static Vector3d? Intersect(Ray3d ray, Plane3d plane, double epsilon = 1e-06)
        {
            // Extract ray and plane components
            var rayOrigin = ray.Origin;
            var rayDirection = ray.Direction.Normalized;
            var planeNormal = plane.Normal.Normalized;
            var planePoint = plane.GetAnyPoint();

            // Calculate the dot product of ray direction and plane normal
            var dotProduct = rayDirection.Dot(planeNormal);

            // Check if the ray is parallel to the plane
            if (Math.Abs(dotProduct) < epsilon)
            {
                // No intersection: ray is parallel to the plane
                return null; 
            }

            // Compute the parameter t for the ray equation
            var t = (planePoint - rayOrigin).Dot(planeNormal) / dotProduct;

            // Check if the intersection lies in the ray's positive direction
            if (t < 0)
            {
                // Intersection is behind the ray origin
                return null; 
            }

            // Compute the intersection point
            var intersectionPoint = rayOrigin + t * rayDirection;
            return intersectionPoint;
        }
    }
}

