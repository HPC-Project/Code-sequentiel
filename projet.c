# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>

//les fonctions de notre programmes

    int main ( int argc, char *argv[] );
    void compute ( int num_particles, int dimension, double position[], double velocity[], double mass, double forces[], double *potential_energy, double *kinetic_energy );
    double cpu_time ();
    double dist ( int dimension, double r1[], double r2[], double diplacement_vector[] );
    void initialize ( int num_particles, int dimension, double position[], double velocity[], double acceleration[] );
    void r8mat_uniform_ab ( int m, int n, double a, double b, int *seed, double r[] );
    void timestamp ();
    void update ( int num_particles, int dimension, double position[], double velocity[], double forces[], double acceleration[], double mass, double size_time_step );

int main ( int argc, char *argv[] )
{
  double *acceleration;
  double ctime;
  double size_time_step;
  double e0;
  double *force;
  double kinetic;
  double mass = 1.0;
  int dimension;
  int num_particles;
  double *position;
  double potential;
  int step;
  int step_num;
  int step_print;
  int step_print_index;
  int step_print_num;
  double *velocity;

  //Important pour comprendre le code: 

/*
Une simulation de dynamique moléculaire consiste à calculer l'évolution d'un système de particules au cours du temps.
 Ces simulations servent de modèles structuraux et dynamiques pour la compréhension de résultats expérimentaux.
*/


/*
Dans ces simulations, le temps évolue de manière discrète:  : le temps est découpé en une suite d'instants  séparés par un intervalle très court appelé "pas-de-temps" ou "time-step"
La simulation consistera alors à calculer la position et la vitesse des particules à chacun des instants, en utilisant les résultats obtenus à l'instant précédent.

Le calcul des forces d'interaction entre les particules permet de déterminer l'évolution des vitesses, et donc des positions, en utilisant les lois de la dynamique classique de Newton discrétisées. 
*/


//d'abord pour pouvoir faire la simulation, on doit avoir la dimension, le nombre de particules et le nombre de time steps

/*
  Get the spatial dimension.
*/
  if ( 1 < argc )
  {
    dimension = atoi ( argv[1] );
  }
  else
  {
    printf ( "\n" );
    printf ( "Entrer la dimension spatiale (2 ou 3).\n" );
    scanf ( "%d", &dimension );
  }
//
//  Get the number of particles.
//
  if ( 2 < argc )
  {
    num_particles = atoi ( argv[2] );
  }
  else
  {
    printf ( "\n" );
    printf ( "Entrer le nombre des particules (500 par exemple).\n" );
    scanf ( "%d", &num_particles );
  }
//
//  Get the number of time steps.
//
  if ( 3 < argc )
  {
    step_num = atoi ( argv[3] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Entrer le nombre du time steps (500 ou 1000 par exemple).\n" );
    scanf ( "%d", &step_num );
  }
//
//  Get the time steps.
//
  if ( 4 < argc )
  {
    size_time_step = atof ( argv[4] );
  }
  else
  {
    printf ( "\n" );
    printf ( "  Entrer le size du time_step (0.1 par exemple).\n" );
    scanf ( "%lf", &size_time_step );
  }


/*
  Allocate memory.
*/

  acceleration = ( double * ) malloc ( dimension * num_particles * sizeof ( double ) );
  force = ( double * ) malloc ( dimension * num_particles * sizeof ( double ) );
  position = ( double * ) malloc ( dimension * num_particles * sizeof ( double ) );
  velocity = ( double * ) malloc ( dimension * num_particles * sizeof ( double ) );

/*
  This is the main time stepping loop:
    Compute forces and energies,
    Update positions, velocities, accelerations.
*/
  printf ( "\n" );
  printf ( "  Dans chaque step, on calcule les energies potentielles et de kinetic.\n" );
  printf ( "  La somme de ces energies doit etre une constante.\n" );
  printf ( "  On affiche aussi l'erreur relative pour vérifier l'accuracy\n" );
  printf ( "  Energie totale\n" );
  printf ( "\n" );
  printf ( "      Step      Energie              Energie        (P+K-E0)/E0\n" );
  printf ( "                Potentielle P        Kinetic K       L'erreur de l'energie relative\n" );
  printf ( "\n" );

  step_print = 0;
  step_print_index = 0;
  step_print_num = 10;

  ctime = cpu_time ( );

//Pour chaque time step: 

  for ( step = 0; step <= step_num; step++ )
  {
      //si step == 0 on doit initialiser la position, la velocité (la vitesse) et l'acceleration de chaque particule --> check la fonction initialize en bas

    if ( step == 0 )
    {
      initialize ( num_particles, dimension, position, velocity, acceleration );
    }
      // sinon: on update les patricules avec les nouvelles positions, vitesses et accelerations --> check la fonction update en bas
    else
    {
      update ( num_particles, dimension, position, velocity, force, acceleration, mass, size_time_step );
    }
     //on calcule les energies et les forces d'interaction entre les particules --> check la fonction compute en bas 
    compute ( num_particles, dimension, position, velocity, mass, force, &potential, &kinetic );

    if ( step == 0 )
    {
      e0 = potential + kinetic;
    }

    if ( step == step_print )
    {
      printf ( "  %8d  %14f  %14f  %14e\n", step, potential, kinetic,
       ( potential + kinetic - e0 ) / e0 );
      step_print_index = step_print_index + 1;
      step_print = ( step_print_index * step_num ) / step_print_num;
    }

  }

//Le temps d'éxecusion total:
  ctime = cpu_time ( ) - ctime;
  printf ( "\n" );
  printf ( "  Le temps d'execution: %f seconds.\n", ctime );
/*
  Free memory.
*/
  free ( acceleration );
  free ( force );
  free ( position );
  free ( velocity );
  return 0;
}

//Le role de cette fonction est d'initialiser la position, la vitesse et l'acceleration de chaque particule 

/* Parameters:

    Input:

    int num_particles: le nombre des particules.

    int dimension: la dimension spatiale.

    Output: 

    double position[dimension*num_particles]: les positions.

    double velocity[dimension*num_particles]: les velocities (les vitesses)

    double acceleration[dimension*num_particles]: les accelerations.

*/
void initialize ( int num_particles, int dimension, double position[], double velocity[], double acceleration[] )
{
  int i;
  int j;
  int seed;
/*
  Set positions.
*/
  seed = 123456789;
  r8mat_uniform_ab ( dimension, num_particles, 0.0, 10.0, &seed, position );
/*
  Set velocities.
*/
  for ( j = 0; j < num_particles; j++ )
  {
    for ( i = 0; i < dimension; i++ )
    {
      velocity[i+j*dimension] = 0.0;
    }
  }
/*
  Set accelerations.
*/
  for ( j = 0; j < num_particles; j++ )
  {
    for ( i = 0; i < dimension; i++ )
    {
      acceleration[i+j*dimension] = 0.0;
    }
  }
  return;
}

//Le role de cette fonction est de mettre à jour la position, la vitesse et l'acceleration de chaque particule 

/* Parameters:

    Input:

    int num_particles: le nombre des particules.

    int dimension: la dimension spatiale.

    double forces[dimension*num_particles]: les forces.

    double mass: la masse.

    double size_time_step: le size du time step.

    Input/output: 
    
    double position[dimension*num_particles]: les positions.

    double velocity[dimension*num_particles]: les velocities (les vitesses)

    double acceleration[dimension*num_particles]: les accelerations.

*/
void update ( int num_particles, int dimension, double position[], double velocity[], double forces[],
  double acceleration[], double mass, double size_time_step )

{
  int i;
  int j;
  double rmass;

  rmass = 1.0 / mass;

  for ( j = 0; j < num_particles; j++ )
  {
    for ( i = 0; i < dimension; i++ )
    {
      position[i+j*dimension] = position[i+j*dimension] + velocity[i+j*dimension] * size_time_step + 0.5 * acceleration[i+j*dimension] * size_time_step * size_time_step;
      velocity[i+j*dimension] = velocity[i+j*dimension] + 0.5 * size_time_step * ( forces[i+j*dimension] * rmass + acceleration[i+j*dimension] );
      acceleration[i+j*dimension] = forces[i+j*dimension] * rmass;
    }
  }

  return;
}

//Le role de cette fonction est de calculer les forces et les energies de chaque patricule
/*  Discussion:

    The computation of forces and energies is fully parallel.

    The potential function V(X) is a harmonic well which smoothly
    saturates to a maximum value at PI/2:

      v(x) = ( sin ( min ( x, PI/2 ) ) )^2

    The derivative of the potential is:

      dv(x) = 2.0 * sin ( min ( x, PI/2 ) ) * cos ( min ( x, PI/2 ) )
            = sin ( 2.0 * min ( x, PI/2 ) )

*/

/*

  Parameters:

    Input: 
    
    int num_particles: le nombre des particules.

    int dimension: la dimension spatiale.

    double position[dimension*num_particles]: les positions.

    double velocity[dimension*num_particles]: les velocities (les vitesses)

    double mass: la masse.

    Output:
    
    double forces[dimension*num_particles]: les forces.

    double *potential_energy: l'energie potentielle totale.

    double *kinetic_energy: l'energie kinetic totale.

*/

void compute ( int num_particles, int dimension, double position[], double velocity[], double mass,
  double forces[], double *potential_energy, double *kinetic_energy )
{
  double d;
  double d2;
  int i;
  int j;
  int k;
  double ke;
  double pe;
  double PI2 = 3.141592653589793 / 2.0;
  double rij[3];

  pe = 0.0;
  ke = 0.0;

  for ( k = 0; k < num_particles; k++ )
  {
/*
  Compute the potential energy and forces.
*/
    for ( i = 0; i < dimension; i++ )
    {
      forces[i+k*dimension] = 0.0;
    }

    for ( j = 0; j < num_particles; j++ )
    {
      if ( k != j )
      {
        d = dist ( dimension, position+k*dimension, position+j*dimension, rij );
/*
  Attribute half of the potential energy to particle J.
*/
        if ( d < PI2 )
        {
          d2 = d;
        }
        else
        {
          d2 = PI2;
        }

        pe = pe + 0.5 * pow ( sin ( d2 ), 2 );

        for ( i = 0; i < dimension; i++ )
        {
          forces[i+k*dimension] = forces[i+k*dimension] - rij[i] * sin ( 2.0 * d2 ) / d;
        }
      }
    }
/*
  Compute the kinetic energy.
*/
    for ( i = 0; i < dimension; i++ )
    {
      ke = ke + velocity[i+k*dimension] * velocity[i+k*dimension];
    }
  }

  ke = ke * 0.5 * mass;

  *potential_energy = pe;
  *kinetic_energy = ke;

  return;
}

//Le role de cette fonction est de calculer le temps d'éxecusion d'un programme
double cpu_time ( )
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}

//Le role de cette fonction est de calculer la distance (et sa norme) entre deux particules
/*
Parameters:

    Input: 
    
    int dimension: la dimension spatiale.

    double R1[dimension], R2[dimension]: les positions des particles.

    Output
    
    double diplacement_vector[dimension]: le vecteur du deplacement.

    double norme_dipalcement, la norme euclidienne du deplacement

*/
double dist ( int dimension, double r1[], double r2[], double diplacement_vector[] )

{
  double norme_dipalcement;
  int i;

  norme_dipalcement = 0.0;
  for ( i = 0; i < dimension; i++ )
  {
    diplacement_vector[i] = r1[i] - r2[i];
    norme_dipalcement = norme_dipalcement + diplacement_vector[i] * diplacement_vector[i];
  }
  norme_dipalcement = sqrt ( norme_dipalcement );

  return norme_dipalcement;
}

//Le role de cette fonction est de retourner a scaled pseudorandom R8MAT.
/* Cette fonction implémente la recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits, including a sign bit.
void r8mat_uniform_ab ( int m, int n, double a, double b, int *seed, double r[] )
*/

/*
 Parameters:

    Input: 
    
    int M, N, the number of rows and columns.

    double A, B, the limits of the pseudorandom values.

    Input/output: 
    int *SEED, the "seed" value.  Normally, this value should not be 0.  On output, SEED has been updated.

    Output: 
    
    double R[M*N], a matrix of pseudorandom values.

*/

void r8mat_uniform_ab ( int m, int n, double a, double b, int *seed, double r[] )

{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
      r[i+j*m] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return;
}

void timestamp ( )
//Le role de cette fonction est d'afficher le temps actuel as a time step
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
