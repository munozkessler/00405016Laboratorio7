#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "tools.h"
#include "display_tools.h"
#include "sel.h"
#include "assembly.h"

int main(int argc, char *argv[])
{
    char filename[150];
    strcpy(filename,argv[1]);

    vector<Matrix> localKs;
    vector<Vector> localbs;
    Matrix K;
    Vector b;
    Vector T;

    cout << "METODO DE LOS ELEMENTOS FINITOS EN DOS DIMENSIONES CON FUNCIONES DE FORMA LINEALES Y PESOS DE GALERKIN - 0040516\n\n";

    mesh m;
    leerMallayCondiciones(m,filename);
    cout << "Datos obtenidos correctamente\n********************\n";

    crearSistemasLocales(m,localKs,localbs);
    showKs(localKs); showbs(localbs);
    cout << "******************************\n";

    zeroes(K,3*m.getSize(NODES));
    zeroes(b,3*m.getSize(NODES));
    ensamblaje(m,localKs,localbs,K,b);
    showMatrix(K); showVector(b);
    cout << "******************************\n";

    applyDirichlet(m,K,b);
    showMatrix(K); showVector(b);
    cout << "******************************\n";

    zeroes(T,b.size());
    calculate(K,b,T);

    cout << "La respuesta es: \n";
    showVector(T);

    writeResults(m,T,filename);

    return 0;
}