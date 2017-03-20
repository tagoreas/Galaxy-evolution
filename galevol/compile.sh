#!/bin/bash

srclist=($(ls *cpp | grep -v template | grep -v kinv))
objlist=($(ls *cpp | grep -v template | grep -v inv | sed s/'.cpp'/'.o'/g))

rm -f *.o libgalevol.so

for f in "${srclist[@]}"; do
    echo mpic++ -O3 -g -fPIC -Wall -fmessage-length=0 -Wno-unused-variable -c "$f"
         mpic++ -O3 -g -fPIC -Wall -fmessage-length=0 -Wno-unused-variable -c "$f"
done

echo mpic++ -O3 -shared -g -fPIC -o libgalevol.so -Wl,-Bdynamic -g ${objlist[@]} -lopenblas
     mpic++ -O3 -shared -g -fPIC -o libgalevol.so -Wl,-Bdynamic -g ${objlist[@]} -lopenblas
