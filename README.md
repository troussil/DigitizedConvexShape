# DigitizedConvexShape

This is a small [DGtal](https://www.dgtal.org/) demo related to the paper "Geometry of Gauss digitized convex shapes" [[springer](https://link.springer.com/chapter/10.1007/978-3-032-09544-2_2)], [[hal](https://cnrs.hal.science/hal-05355892v1)]. 

# Run 

The program takes a convex shape name and a grid step as input. It computes the Gauss digitization of the shape, extracts the discrete surface and its convex hull, then estimates a normal vector at each surface element using the convex‐hull facet onto which it projects. These estimates are compared with the normals of the input shape. The angular error is displayed as a color map on the discrete surface, and statistics (minimum/average/maximum angular errors) are computed.

For instace, by running `./main -p ellipsoid -g 0.2`, you can see this in the viewer: 

![gif animation showing several snapshots](https://github.com/troussil/DigitizedConvexShape/blob/main/images/ellipsoid.gif)

# Build

Once [DGtal](https://www.dgtal.org/) is built (with the [Polyscope](https://polyscope.run/) option enabled), you can compile as follows: 
- `mkdir build`
- `cd build`
- `cmake .. -DDGtal_dir=<path to your DGtal build directory> -DCMAKE_BUILD_TYPE=Release`
- `make`
