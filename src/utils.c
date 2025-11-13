// gcc -O3 -fPIC -shared -o libutils.dylib utils.c -lm
// gcc -O3 -fPIC -shared -o libutils.so utils.c -lm
// gcc -O3 -shared -o utils.dll utils.c -lm 

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAXLON   -4
#define MAXLAT   70
#define MINLAT   60
#define MINLON  -32
#define PI       3.14159265358979323846264338327950288419716939937510582097494459230781640628

#ifdef _WIN  // Windows platform
  #include <windows.h>
  #define EXPORT __declspec(dllexport)
#else  // POSIX platforms (Linux/macOS)
  #define EXPORT
#endif

double degmin2deg(double degmin) {

  double min ;

  if (fabs(degmin)<10000)
    degmin = degmin*100 ;
  min = (degmin/100)-floor(degmin/10000.)*100. ;
  return(((degmin+(200.0/3.0)*min)/10000)) ;
}

void deg2point (double *x, double *y, double *Lat, double *Lon, int length, int Norm) {

  double lat65, x65, M, M65, Diff, lat, scale, lon ;
  int    i ;

  /* at 65�lat 18�lon */
  lat65 = 65.*PI/180. ;
  x65 = (111415.13*cos(lat65)-94.55*cos(3*lat65)+0.12*cos(5*lat65))/60 ;
  M65 = 7915.704456*log10(tan(PI/4+lat65/2))
           - sin(lat65)*(23.110771+0.052051*(sin(lat65))*sin(lat65)) ;

  /* Automatic Scaling */
  lat = MAXLAT*PI/180 ;
  M = 7915.704456*log10(tan(PI/4+lat/2))-sin(lat)*(23.110771+0.052051*(sin(lat))*sin(lat65)) ;
  Diff = M-M65 ;
  lat = MINLAT*PI/180 ;
  M = 7915.704456*log10(tan(PI/4+lat/2))-sin(lat)*(23.110771+0.052051*(sin(lat))*sin(lat65)) ;
  Diff = Diff + M65 - M ;
  scale = (double)Norm/(Diff*x65) ;
  x65 = scale*x65 ;

  for (i=0; i<length; i++) {
    lon = Lon[i];
    lat = Lat[i]*PI/180 ;
    M = 7915.704456*log10(tan(PI/4+lat/2))-sin(lat)*(23.110771+0.052051*(sin(lat))*sin(lat65)) ;
    Diff = M65-M ;
    y[i] = (int) (Diff*x65) ;
    x[i] = (int) ((lon+18)*x65*60) ;
  }
}

/* Uses the dot product to determine if two line intercept
   I think we should define waypoint along the coast to use for seeing if we cross the land */
int intersects(double s0x[2], double s0y[2], double s1x[2], double s1y[2]) {
  double dx0 = s0x[1] - s0x[0];
  double dx1 = s1x[1] - s1x[0];
  double dy0 = s0y[1] - s0y[0];
  double dy1 = s1y[1] - s1y[0];
  double p0 = dy1 * (s1x[1] - s0x[0]) - dx1 * (s1y[1] - s0y[0]);
  double p1 = dy1 * (s1x[1] - s0x[1]) - dx1 * (s1y[1] - s0y[1]);
  double p2 = dy0 * (s0x[1] - s1x[0]) - dx0 * (s0y[1] - s1y[0]);
  double p3 = dy0 * (s0x[1] - s1x[1]) - dx0 * (s0y[1] - s1y[1]);
  return ((p0 * p1 <= 0) && (p2 * p3 <= 0));
}

int crossesland(double stAx, double  stAy, double stBx, double stBy, double *Land, int n) {
  int i, dn = 4;
  double s0x[2], s0y[2], s1x[2], s1y[2];
  
  /*if ((stAx == stBx) && (stAy == stBy))
    return 1;*/

  s0x[0] = degmin2deg(stAx);
  s0x[1] = degmin2deg(stBx);
  s0y[0] = -degmin2deg(stAy);
  s0y[1] = -degmin2deg(stBy);
  deg2point(s0x,s0y,s0x,s0y,2,1000000);
  
  for (i = 0; i < n-dn; i = i + dn) { 
    s1x[0] = Land[i];
    s1x[1] = Land[i+dn];
    s1y[0] = Land[i+n];
    s1y[1] = Land[i+dn+n];
    deg2point(s1x,s1y,s1x,s1y,2,1000000);
    if (1 == intersects(s0x,s0y,s1x,s1y)) {
      return 0;
    }
  }
  return 1;
}

int createfeasiblelinkmatrix(int *M, double *Data[4], int m, double *LandDeg, int n) {
  int i, j;

  for (i = 0; i < m; i++)
    M[i+m*i] = 1;
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++) {
      M[i+m*j] = crossesland(Data[0][i], Data[1][i], Data[0][j], Data[1][j], LandDeg, n);
      M[j+m*i] = M[i+m*j];
    }
  }
  return 0;
}

int minDistance(double *dist, int *sptSet, int n) {
  double min = 100000000.0; /* Infinity */
  int v, min_index;
 
  for (v = 0; v < n; v++) {
    if ((sptSet[v] == 0) && (dist[v] <= min)) {
      min = dist[v];
      min_index = v;
    }
  }
  return min_index;
}

int dijkstra(double *d, int *path, double **graph, int n, int src, int dest) {
  double *dist, INFTY = 10000000000.0;
  int *sptSet, i, j=0, count, u, v, *nodes;

  dist = (double *) malloc(n*sizeof(double)); 
  sptSet = (int *) calloc(n, sizeof(int));
  nodes = (int *) calloc(n, sizeof(int));
  for (i = 0; i < n; i++)
    dist[i] = INFTY; /* smaller than infinity above */
  dist[src] = 0.0;
  for (count = 0; count < n - 1; count++) {
    u = minDistance(dist, sptSet, n);
    sptSet[u] = 1;
    for (v = 0; v < n; v++) {
      if ((sptSet[v]==0) && (graph[u][v] > 0) && (dist[u] < INFTY) && (dist[u] + graph[u][v] < dist[v])) {
        dist[v] = dist[u] + graph[u][v];
        nodes[v] = u;
      }
    }
  }
  nodes[src] = dest;
  i = dest;
  path[j++] = i;
  while (i != src) {
    i = nodes[i];
    path[j++] = i;
  }
  path[j] = -1; /* end of path */
  *d = dist[dest];
  free(dist);
  free(sptSet);
  free(nodes);
  return 0;
}
