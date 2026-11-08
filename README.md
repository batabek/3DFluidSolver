# 3DFluidSolver
3D Fluid Solver based on Stable Fluids paper.

* This is an implementation of a 3D Fluid Solver based on Jos Stams' "Stable Fluids" paper.
  It uses a colocated grid instead of a staggered grid since this is straightforward to implement.
* Main program outputs the simulation frames in PBRT format.

## Building

```
make        # builds ./FluidSolver (requires a C++17 compiler)
make run    # runs the simulation, frames are written to Out/Frame_N.pbrt
```

## Rendering with pbrt-v4

Each frame is written as a pbrt-v4 `MakeNamedMedium "smoke"` definition. `smoke.pbrt` is a
ready-made scene that pulls a frame in with `Include` and attaches the medium to a unit cube:

```
Include "Out/Frame_100.pbrt"
...
AttributeBegin
    MediumInterface "smoke" ""
    Material "interface"
    Shape "trianglemesh" ...   # unit cube bounding the grid
AttributeEnd
```

pbrt-v4 has no Homebrew package, so build it from source (needs cmake):

```
git clone --recursive https://github.com/mmp/pbrt-v4.git ~/src/pbrt-v4
cmake -S ~/src/pbrt-v4 -B ~/src/pbrt-v4/build -DCMAKE_BUILD_TYPE=Release
cmake --build ~/src/pbrt-v4/build -j --target pbrt_exe
```

Then render:

```
~/src/pbrt-v4/build/pbrt smoke.pbrt                  # frame 100 -> smoke.png
PBRT=~/src/pbrt-v4/build/pbrt ./render.sh            # every frame -> Out/Frame_N.png
PBRT=~/src/pbrt-v4/build/pbrt ./render.sh 10 50      # selected frames
```

The smoke's look (`sigma_a`, `sigma_s`, `scale`) is set in `exportGridToPBRTFormat`.

## Making an animation

After rendering every frame, combine them with ffmpeg (`brew install ffmpeg`):

```
./animate.sh        # Out/Frame_N.png -> smoke.mp4 and smoke.gif at 24 fps
./animate.sh 12     # slower, 12 fps
```
