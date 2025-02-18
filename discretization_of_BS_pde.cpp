// for discretization purposes:
//          crank nicholson
//          finite pieces method
//          the diferential opertaor aplied to the variation V
//  -------- will give a system of tridiagonal eq thus we have to solve it inteligently'
// -------------- Thomas Algo for this purpose

// discretiztion
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// some random options parameters
double  Smax = 100;  // max grile price
double K = 50;     // strike price
double T = 1.0;   // expiration time
double volatility = 0.2;
double stepsInPrice = 100;
double stepsInTime = 1000;

class TridiagonalSolver{
    // so this a class that solves the sysyem of eq with the crank nicholson approach
    public:
        // we solve the tridiagonal sistem of coeff
        static void solve( 
            std::vector<double>& a , std::vector<double>& b ,
            std::vector<double>& c , std::vector<double>& d ,
            std::vector<double>& x
         ){
            int n = b.size();
            std::vector<double> c_star( n , 0 );
            std::vector<double> d_star( n , 0 ); // here the std::vectors behave as constructors 

            // sweeping onward
            c_star[0] = c[0] / b[0];
            d_star[0] = d[0] / b[0];
            for( int i = 1; i < n ; ++ i ){
                double m = 1.0 / (
                    b[i] - a[i] * c_star[i-1]
                );
                c_star[i] = c[i] * m;
                d_star[i] = ( 
                    d[i] - a[i] * d_star[i-1]
                 ) * m;
            }

            //backward substitution
            x[n-1] = d_star[n-1];
            for( int i = n - 2 ; i >= 0 ; i-- ){
                x[i] = d_star[i] - c_star[i] * x[i+1];
            }
         }  
};


// Black Scholes discretization and solving
class BlackScholesSolver{
    private:
        double Smax , K , T , r , sigma;
        int M , N;
        std::vector<std::vector<double>> V; // the volatility is a 2 d tensor a matrix 
    public:
        BlackScholesSolver( double Smax_ , double K_,
            double T_ , double r_ , double sigma_ , int M_ , int N_ ) : 
            Smax( Smax_ ) , K(K_) , T(T_) , r(r_) , sigma(sigma_) , M(M_) , N(N_){
                V.resize( N+1 , vector<double>( M+1 , 0.0 ) );
                // re-arranging the matrix
            }
        void initialize(){
            double dS = Smax / M;
            for( int i = 0 ; i <= M ; i++ ){
                double S = i * dS;
                V[N][i] = max( S - K , 0.0 ); // expiration payoff
            }
        }
        void solve(){
            double dS = Smax / M;
            double dt = T / N;
            std::vector<double> a(M-1,0), b(M-1,0), c(M-1,0), d(M-1,0) , x(M-1,0);
            for( int n = N-1 ; n >= 0 ; n-- ){
                for( int i = 1 ; i < M ; ++ i){
                    double S = i * dS;
                    double alphaCoeff = 0.25 * dt * ( pow(sigma,2) * pow(S,2) / pow( dS , 2 ) - r * S /dS );
                    double betaCoeff = -0.5 * dt * ( pow(sigma,2) * pow(S,2) / pow(dS,2) + r );
                    double gammaCoeff = 0.25 * dt * ( pow(sigma,2) * pow(S,2) / pow(dS,2) + r*S/dS );
                
                    if( i>1 ) a[i-1] = -alphaCoeff;
                    b[i-1] = 1 - betaCoeff;
                    if( i<M-1 ) c[i-1] = -gammaCoeff;
                    d[i-1] = alphaCoeff*V[n+1][i-1] + 
                        (1+betaCoeff)*V[n+1][i] + gammaCoeff * V[n+1][i+1];
                }
                TridiagonalSolver::solve( a,b,c,d,x );
                for( int i=1; i<M ; ++i ){
                    V[n][i] = x[i-1];
                }
            }
        }
        double get_option_price( double S0 ){
            int index = static_cast<int>(S0/(Smax/M));
            return V[0][index];
        }
};


int main(){
    double SMax = 100 , K = 50 , T = 1.0 , r = 0.05 , sigma = 0.2;
    int M = 100 , N = 1000;
    BlackScholesSolver  BSs( SMax , K , T , r , sigma , M , N );
    BSs . initialize();
    BSs.solve();

    double S0 = 50;
    std::cout << "estimated value for option S0 = " << S0 << " is ---- " << BSs .get_option_price( S0 ) << '\n';

return 0;
}