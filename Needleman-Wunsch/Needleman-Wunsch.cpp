/*
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
*/

#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <map>
#include <vector>

static int * SymbolToIndex = NULL ;
static int * Blosum62 = NULL ;

static void LoadBlosum62( const char * Path )
{
    std::map< char , int > Symbols ;
    std::vector< int > Values ;

    std::ifstream File( Path ) ;
    std::string Line ;


    int Index = 0 ;
    while ( std::getline( File , Line ) )
    {
        if ( Line[0] == '#' )
        {
            continue ;
        }

        std::istringstream Input( Line ) ;
        char Symbol = '*' ;
        Input >> Symbol ;
        if ( Symbol == '*' )
        {
            Symbol = '[' ;
        }
        Symbols.insert( std::make_pair( Symbol , Index ++ ) ) ;

        int Value = - 4 ;
        while ( Input >> Value )
        {
            Values.push_back( Value ) ;
        }
    }

    SymbolToIndex = new int[26] ;
    std::fill( SymbolToIndex , SymbolToIndex + 26 , - 1 ) ;
    for ( std::map< char , int >::iterator itr = Symbols.begin() ; itr != Symbols.end() ; ++ itr )
    {
        SymbolToIndex[itr->first - 'A'] = itr->second ;
    }

    Blosum62 = new int[Values.size()] ;
    std::copy( Values.begin() , Values.end() , Blosum62 ) ;
}

int LookupBlosum( char A , char B )
{
    int X = SymbolToIndex[A - 'A'] ;
    int Y = SymbolToIndex[B - 'A'] ;

    int V = Blosum62[ Y * 24 + X ] ;

    return V ;
}

int main( int Argc , char * Argv[] )
{
    -- Argc , ++ Argv ;
    if ( Argc != 3 )
    {
        return 1 ;
    }

    //
    LoadBlosum62( Argv[0] ) ;

    //
    char * A = Argv[1] ;
    char * B = Argv[2] ;

    //
    int W = strlen( A ) + 1 ;
    int H = strlen( B ) + 1 ;

    int * Score = new int[H * W] ;
    std::fill( Score , Score + H * W , 0 ) ;

    char * Trace = new char[H * W] ;
    std::fill( Trace , Trace + H * W , 0 ) ;

    for ( int x = 0 ; x < W ; ++ x )
    {
        Score[x] = - x * 10 ;
        Trace[x] = 'L' ;
    }
    for ( int y = 0 ; y < H ; ++ y )
    {
        Score[y * W] = - y * 10 ;
        Trace[y * W] = 'U' ;
    }
    Trace[0] = 'E' ;

    //
    for ( int y = 1 ; y < H ; ++ y )
    {
        for ( int x = 1 ; x < W ; ++ x )
        {
            int S = LookupBlosum( A[x - 1] , B[y - 1] ) ;
            int D = Score[( y - 1 ) * W + ( x - 1 )] + S ;
            int U = Score[( y - 1 ) * W + x] - 10 ;
            int L = Score[y * W + x - 1] - 10 ;

            int V = 0 ;
            char T = 'X' ;
            if ( D >= U && D >= L )
            {
                V = D ;
                T = 'D' ;
            }
            else if ( U >= D && U >= L )
            {
                V = U ;
                T = 'U' ;
            }
            else if ( L >= D && L >= U )
            {
                V = L ;
                T = 'L' ;
            }
            else
            {
                printf( "(%d %d)[%d %d %d]\n" , x , y , D , U , L ) ;
            }

            Score[y * W + x] = V ;
            Trace[y * W + x] = T ;
        }
    }

    //

#if 0
    for ( int y = 0 ; y < H ; ++ y )
    {
        for ( int x = 0 ; x < W ; ++ x )
        {
            printf( "%d:%c " , Score[y * W + x] , Trace[y * W + x] ) ;
        }
        std::cout << std::endl ;
    }
#endif

    //
    std::vector< char > Aligned ;
    int x = W - 1 , y = H - 1 ;
    while ( (x > 0) && (y > 0) )
    {
        char V = Trace[y * W + x] ;
        switch ( V )
        {
            case 'D' :
                Aligned.push_back( B[y - 1] ) ;
                -- x , -- y ;
                break ;
            case 'U' :
                Aligned.push_back( '-' ) ;
                -- y ;
                break ;
            case 'L' :
                Aligned.push_back( '-' ) ;
                -- x ;
                break ;
            default :
                break ;
        }
    }
    std::cout << "Result: " << std::endl ;
    std::copy( Aligned.rbegin() , Aligned.rend() , std::ostream_iterator< char >( std::cout , "" ) ) ;
    std::cout << std::endl ;


    delete [] Score ;
    delete [] Trace ;

	return 0 ;
}