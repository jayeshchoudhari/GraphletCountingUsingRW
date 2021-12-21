# Counting graphlets using Random walk

## Input file format:
```
#nodes #edges
v1 v2
v2 v4
v3 v4
.
.
.
```
## Graphlet codes:
1. **g32**: Triangle (3-Clique)\n
2. **g43**: 4-Cycle\n
3. **g45**: 4-Chord-Cycle\n
4. **g46**: 4-Clique\n
5. **g58**: 5-Clique-But-2-Edges-(Hatted-4-Clique)\n
6. **g59**: 5-Clique-But-1-Edge-(Almost-5-Clique)\n
7. **g510**: 5-Clique\n
8. **g615**: 6-Clique\n

## Compile the code:
```
make
``` 

### Run the code:
Running the code requires 4 parameters:
1. Input Graph filename 
2. CodeWord for Graphlet To Count 
3. Graph name and
4. Random walk edges(0) or U.A.R. edges (1)

### To count 3-Clique or Triangles using random walk edges
```
./main-DiffStartPoint orkut.edges g32 orkut 0
```