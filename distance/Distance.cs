﻿namespace g4
{
    // [TODO] get rid of this!!
    public static class Distance
    {
        public static float ClosestPointOnLineT(Vector3f p0, Vector3f dir, Vector3f pt)
        {
            float t = (pt - p0).Dot(dir);
            return t;
        }
    }
}
