# include "tsunami.h"

//definition des structure
typedef struct {
    int elem[2];
    int node[2];
    int nEdge;
    int number;
} femEdge;

typedef struct {
    int* elem;
    double* X;
    double* Y;
    double* bath;
    int nElem;
    int nNode;
} femMesh;

typedef struct {
    femMesh* mesh;
    femEdge* edges;
    int size;
    double* E;
    double* U;
    double* V;
    double* FE;
    double* FU;
    double* FV;
} femTsunamiProblem;

//definition des fonction
double interpolate(double* phi, double* U, int* map, int n);
void femTsunamiAddIntegralsElements(femTsunamiProblem* myTsunami);

//fonction principale
void tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName)
{
    femMesh* theMesh = malloc(sizeof(femMesh));
    femEdge* theEdges = malloc(sizeof(femEdge));
    int i, trash, * elem;

//lecture du Fichier
    FILE* file = fopen(meshFileName, "r");             
    if (file == NULL) printf("No mesh file !");
    fscanf(file, "Number of nodes %d \n", &theMesh->nNode);
    theMesh->X = malloc(sizeof(double) * theMesh->nNode);
    theMesh->Y = malloc(sizeof(double) * theMesh->nNode);
    theMesh->bath = malloc(sizeof(double) * theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        fscanf(file, "%d : %le %le %le \n", &trash, &theMesh->X[i], &theMesh->Y[i]), & theMesh->bath;
    }
    fscanf(file, "Number of triangles % d \n", &theMesh->nElem);
    theMesh->elem = malloc(sizeof(int) * 3 * theMesh->nElem);
    for (i = 0; i < theMesh->nElem; ++i) {
        elem = &(theMesh->elem[i * 3]);
        fscanf(file, "%d : %d %d %d\n", &trash, &elem[0], &elem[1], &elem[2]);
    }
    fscanf(file, "Number of edges %d \n", &theEdges->nEdge);
    for (i = 0; i < theEdges->nEdge; i++) {
        fscanf(file, "%d : %d %d %d %d \n", &trash, &theEdges[i].node[0], &theEdges[i].node[1], &theEdges[i].elem[0], &theEdges[i].elem[1]);
    }
    fclose(file);

    femTsunamiProblem* myTsunami = malloc(sizeof(femTsunamiProblem)); //alocation de mémoire pour mon probleme
    myTsunami->mesh = theMesh;
    myTsunami->edges = theEdges;
    int size = myTsunami->mesh->nElem * 3 + 1;
    myTsunami->size = size;
    myTsunami->E = malloc(sizeof(double) * size);
    myTsunami->U = malloc(sizeof(double) * size);
    myTsunami->V = malloc(sizeof(double) * size);
    myTsunami->FE = malloc(sizeof(double) * size);
    myTsunami->FU = malloc(sizeof(double) * size);
    myTsunami->FV = malloc(sizeof(double) * size); 
    for (int i = 0; i < myTsunami->size; i++) {
        myTsunami->FE[i] = 0.0;
        myTsunami->FU[i] = 0.0;
        myTsunami->FV[i] = 0.0;
    }
    femMesh* mesh = myTsunami->mesh;
     
    //Condition initial de Okada
    for (i = 0; i < mesh->nElem; ++i) {
        int* Coord = &(mesh->elem[3 * i]);
        for (int j = 0; j < 3; ++j) {
            myTsunami->E[3 * i + j] = tsunamiInitialConditionOkada(mesh->X[Coord[j]], mesh->Y[Coord[j]]);
            myTsunami->U[3 * i + j] = 0.0;
            myTsunami->V[3 * i + j] = 0.0;
        }
    }




    /*int j, nElem;
    double dtrash;
    
    double BathMax = 9368;
    fscanf(file, "Number of nodes %d \n",&nNode);   
    double *bath = malloc(sizeof(double)*nNode);
    for (i = 0; i < nNode; i++) 
        fscanf(file,"%d : %le %le %le\n",&trash,&dtrash,&dtrash,&bath[i]); 
    fscanf(file, "Number of triangles %d \n",&nElem); 
    int *elem = malloc(sizeof(int)*3*nElem);
    for (i = 0; i < nElem; i++) 
        fscanf(file,"%d : %d %d %d \n",&trash,&elem[i*3],&elem[i*3+1],&elem[i*3+2]);   
    fclose(file); 
    
    double *E  = malloc(sizeof(double)*nElem*3);
    for (i = 0; i < nElem; i++)
        for (j = 0; j < 3; j++)
            E[i*3+j] = bath[elem[i*3+j]]/(10*BathMax);
  
    tsunamiWriteFile(baseResultName,0,E,E,E,nElem,3); 
    
    free(bath);
    free(E);
    free(elem);*/
 
}

void femTsunamiAddIntegralsElements(femTsunamiProblem* myTsunami)
{
    double  xLoc[3], yLoc[3], phi[3], dphidx[3], dphidy[3];
    double  xsi, eta, weight, jac;
    double  y, x, bath, e, u, v;
    int     i, j, k, mapElem[3];
    double* BE = myTsunami->FE;
    double* BU = myTsunami->FU;
    double* BV = myTsunami->FV;
    double* E = myTsunami->E;
    double* U = myTsunami->U;
    double* V = myTsunami->V;
    double* Y = myTsunami->mesh->Y;
    double* X = myTsunami->mesh->X;
    double* Bath = myTsunami->mesh->bath;
    int* nElem = myTsunami->mesh->nElem;
    int* elem = myTsunami->mesh->elem;

    for (i = 0; i < nElem; i++) {
       
        int* mapCoord = &(elem[i * 3]);
        for (j = 0; j < 3; ++j) {
            xLoc[j] = X[mapCoord[j]];
            yLoc[j] = Y[mapCoord[j]];
        }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        dphidx[0] = (yLoc[1] - yLoc[2]) / jac;
        dphidx[1] = (yLoc[2] - yLoc[0]) / jac;
        dphidx[2] = (yLoc[0] - yLoc[1]) / jac;
        dphidy[0] = (xLoc[2] - xLoc[1]) / jac;
        dphidy[1] = (xLoc[0] - xLoc[2]) / jac;
        dphidy[2] = (xLoc[1] - xLoc[0]) / jac;

        for (j = 0; j < 3; j++) {

            xsi = gaussTriangleXsi[j];
            eta = gaussTriangleXsi[j];
            weight = gaussTriangleWeight[j];
            y = interpolate(phi, Y, mapCoord, 3);
            x = interpolate(phi, X, mapCoord, 3);
            bath = interpolate(phi, Bath, mapCoord, 3);
            e = interpolate(phi, E, mapElem, 3);
            u = interpolate(phi, U, mapElem, 3);
            v = interpolate(phi, V, mapElem, 3);
            
            const double Rapport_longeur = ((4 * R * R + x * x + y * y) / (4 * R * R));
            const double f = 2 * Omega * ((4 * R * R - x * x - y * y) / (4 * R * R));

            for (i = 0; i < 3; i++)
            {
                BE[mapElem[i]] += ((dphidx[i] * bath * u + dphidy[i] * bath * v) * Rapport_longeur + (phi[i] * (bath * (x * u + y * v) / (R * R)))) * jac * weight;
                BU[mapElem[i]] += ((phi[i] * (f * v - Gamma * u) + dphidx[i] * g * e * Rapport_longeur) + (phi[i] * g * x * e / (2 * R * R))) * jac * weight;
                BV[mapElem[i]] += ((phi[i] * (-f * u - Gamma * v) + dphidy[i] * g * e * Rapport_longeur) + (phi[i] * g * y * e / (2 * R * R))) * jac * weight;
            }
        }
    }
}

void femTsunamiAddIntegralsEdges(femTsunamiProblem* myTsunami)  //pas fini
{
    double* BE = myTsunami->FE;
    double* BU = myTsunami->FU;
    double* BV = myTsunami->FV;
    double* E = myTsunami->E;
    double* U = myTsunami->U;
    double* V = myTsunami->V;
    femEdge* edge = myTsunami->edges;
    femMesh* mesh = myTsunami->mesh;

    double  xEdge[2], yEdge[2], phiEdge[2], bathEdge[2];
    double  xsi, weight, jac;
    double  eL, eR, uL, uR, vL, vR, unL, unR;
    double  qe, qu, qv;
    int     i, j,k, mapEdge[2][2];
    int     sizeGlo = mesh->nElem * 3 + 1;

    for (i = 0; i < edge->nEdge; i++) {
        for (j = 0; j < 2; ++j) {
            int node = edge[i].node[j];
            for (k = 0; k < 2; k++) {
                int elem = edge[i].elem[k];
                mapEdge[k][j] = (mesh->nElem) * 3;
                if (elem >= 0) {
                    for (i = 0; i < 3; i++) {
                        if (mesh->elem[elem * 3 + i] == node) {
                            mapEdge[k][j] = elem * 3 + i;
                        }
                    }
                }
            }
        }
        for (j = 0; j < 2; ++j) {
            int node = edge[i].node[j];
            xEdge[j] = mesh->X[i];
            yEdge[j] = mesh->Y[i];
            bathEdge[j] = mesh->bath[node];
        }

        int boundary = (mapEdge[1][0] == sizeGlo - 1);

        double dxdxsi = (xEdge[1] - xEdge[0]);
        double dydxsi = (yEdge[1] - yEdge[0]);
        double norm = sqrt(dxdxsi * dxdxsi + dydxsi * dydxsi);
        double nx = dydxsi / norm;
        double ny = -dxdxsi / norm;
        jac = norm / 2.0;
        for (k = 0; k < theRule->n; k++) {
            xsi = theRule->xsi[k];
            weight = theRule->weight[k];
            femDiscretePhi1(theSpace, xsi, phiEdge);
            eL = interpolate(phiEdge, E, mapEdge[0], 2);
            eR = boundary ? eL : interpolate(phiEdge, E, mapEdge[1], 2);
            uL = interpolate(phiEdge, U, mapEdge[0], 2);
            uR = interpolate(phiEdge, U, mapEdge[1], 2);
            vL = interpolate(phiEdge, V, mapEdge[0], 2);
            vR = interpolate(phiEdge, V, mapEdge[1], 2);
            unL = uL * nx + vL * ny;
            unR = boundary ? -unL : uR * nx + vR * ny;
            qe = 0.5 * h * ((unL + unR) + sqrt(g / h) * (eL - eR));
            qu = 0.5 * g * nx * ((eL + eR) + sqrt(h / g) * (unL - unR));
            qv = 0.5 * g * ny * ((eL + eR) + sqrt(h / g) * (unL - unR));
            for (i = 0; i < 2; i++) {
                BE[mapEdge[0][i]] -= qe * phiEdge[i] * jac * weight;
                BU[mapEdge[0][i]] -= qu * phiEdge[i] * jac * weight;
                BV[mapEdge[0][i]] -= qv * phiEdge[i] * jac * weight;
                BE[mapEdge[1][i]] += qe * phiEdge[i] * jac * weight;
                BU[mapEdge[1][i]] += qu * phiEdge[i] * jac * weight;
                BV[mapEdge[1][i]] += qv * phiEdge[i] * jac * weight;
            }
        }
    }
}

double interpolate(double* phi, double* U, int* map, int n) {
    double u = 0.0; int i;
    for (i = 0; i < n; i++)
        u += phi[i] * U[map[i]];
    return u;
}













