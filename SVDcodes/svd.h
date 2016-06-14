
int readNSymA(char *filename, double ***A, int *cols);
double calcNorm(double **A, int rows, int cols);
double calcOffA(double **A, int rows, int cols);
void rowRot(double ***A, int cols, int p, int q, double c, double s);
void colRot(double ***A, int rows, int p, int q, double c, double s);
void JacobiCS(double Apq, double App, double Aqq, double &c, double &s, double tol);
void RandJacobiCS(double Apq, double App, double Aqq, double &c, double &s, double x, double tol);
unsigned long CyclicJacobi(double ***A, int n, double eps, double tol, double param);
unsigned long IndependentJacobi(double ***A, int n, double eps, double tol);
unsigned long RandomJacobi(double ***A, int n, double eps,  double tol, double randParam);
unsigned long BlockRandomJacobi(double ***A, int n, double eps, double tol, double randParam);
void music(int *top, int *bot, int **newtop, int **newbot, int m);
unsigned long CyclicOneJacobi(double ***A, int m, int n, double eps, double tol, double param);
unsigned long IndependentOneJacobi(double ***A, int m, int n, double eps, double tol);
unsigned long RandomOneJacobi(double ***A, int m, int n, double eps, double tol, double randParam);
unsigned long SortedOneJacobi(double ***A, int m, int n, double eps, double tol, double randParam);
unsigned long BlockRandomOneJacobi(double ***A, int m, int n, double eps, double tol, double randParam);
unsigned long StrumpenJacobi(double ***A, int m, int n, double eps, double tol, double randParam, int R);
unsigned long StrumpenRelaxationJacobi(double ***A, int m, int n, double eps, double tol, double randParam, int R);

double vectornorm(double **A, int p, int n);
double dotproduct(double **A, int p, int q, int n);
