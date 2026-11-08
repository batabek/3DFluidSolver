#!/bin/sh
# Combines the rendered frames Out/Frame_N.png into smoke.mp4 and smoke.gif with ffmpeg.
# usage: ./animate.sh [fps]   (default 24)
set -e
cd "$(dirname "$0")"

FPS=${1:-24}

if ! command -v ffmpeg >/dev/null 2>&1; then
    echo "ffmpeg is not installed, install it with: brew install ffmpeg" >&2
    exit 1
fi

# frames are numbered from 1 without zero padding, which %d matches
ffmpeg -y -loglevel error -framerate "$FPS" -start_number 1 -i Out/Frame_%d.png \
    -c:v libx264 -pix_fmt yuv420p -crf 18 smoke.mp4
echo "Wrote smoke.mp4"

# build a palette from the frames first so the gif keeps smooth gradients
ffmpeg -y -loglevel error -framerate "$FPS" -start_number 1 -i Out/Frame_%d.png \
    -vf "split[a][b];[a]palettegen[p];[b][p]paletteuse" smoke.gif
echo "Wrote smoke.gif"
