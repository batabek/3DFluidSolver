#!/bin/sh
# Renders simulated frames with pbrt-v4 using smoke.pbrt as the scene.
# usage: ./render.sh            render every frame in Out/
#        ./render.sh 10 50 100  render only the given frames
# Set PBRT to the pbrt executable if it is not on your PATH.
set -e
cd "$(dirname "$0")"

PBRT=${PBRT:-pbrt}

if [ $# -eq 0 ]; then
    set -- $(ls Out/Frame_*.pbrt | sed 's|Out/Frame_||; s|\.pbrt||' | sort -n)
fi

for i in "$@"; do
    # scene files resolve Include paths relative to themselves, so put each one next to its frame
    sed "s|\"Out/Frame_100.pbrt\"|\"Frame_$i.pbrt\"|" smoke.pbrt > "Out/Scene_$i.pbrt"
    "$PBRT" --quiet --outfile "Out/Frame_$i.png" "Out/Scene_$i.pbrt"
    echo "Rendered Out/Frame_$i.png"
done
