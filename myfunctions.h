void Print_Matrix(double A[], int n, int m);

void Print_Vector(double b[], int n);

void Build_K_elemental(double Ke[], int n, double A, double E, double l, double I);

void Build_F_elemental(double Fe[], double qx, double qy, double l, double time);

void Build_global_matrix(double K[], double Ke[], int Nx, int N);

void Build_K_global_banded(double Kb[], double Ke[], int Nx, int N, int kl, int ku);

void Banded_Storage(double Kb[], double K[], int N, int kl, int ku);

void Build_F_global(double F[], double Fe[], int Nx, double Fy, double time);

void Matrix_System_Solver(int N, double K[], double F[]);

void Write_Text_File(double F[], int N, double l, double L);

void Build_M_elemental(double Me[], int Ne, double rho, double A, double l);

void Build_Un1_Multiplier(double K_temp[], double K[], double M[], int N, double del_t);

void Build_Multiplier(double S[], double F[], double K_temp[], double M[], double del_t, double u0[], double u1[], int N);

void Copy_Vector(int N, double a[], double b[]);