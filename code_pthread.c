# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <pthread.h>

# define NUM_THREADS 5
pthread_barrier_t barrier;
pthread_mutex_t mutex_step;
pthread_mutex_t mutex_step_print;
pthread_mutex_t mutex_ke;
pthread_mutex_t mutex_pe;
pthread_mutex_t mutex_pe_ke;
int dimension;
int num_particles;
double *position;
double *acceleration;
double *velocity;
double *forces;
double e0;
double kinetic;
double mass = 1.0;
double potential;
double size_time_step;
int step; //l'indice d'xecution des seteps
int step_num; 
int step_print;
int step_print_index; //le nombre des prints
int step_print_num; //le pas d'affichage

//les fonctions de notre programmes

    int main ( int argc, char *argv[] );
    void* parallel(void* tid);
    void *compute (int tid);
    double cpu_time ();
    double dist (double r1[], double r2[], double diplacement_vector[]);
    void initialize ();
    void r8mat_uniform_ab ( int m, int n, double a, double b, int *seed, double r[] );
    void timestamp ();
    void *update (int ii);

//   typedef struct initialize_params{
//   long* i; 
//   int num_particles; 
//   int dimension; 
//   double *position; 
//   double *velocity; 
//   double *acceleration;
// }initialize_params;

// typedef struct compute_params{
//   long* i;
//   int* num_particles; 
//   int dimension;
//   double *position;
//   double *velocity;
//   double mass;
//   double *forces;
//   double *potential_energy; 
//   double *kinetic_energy;
// } compute_params;

// typedef struct update_params{
//   long* i;
//   int num_particles; 
//   int dimension;
//   double *position;
//   double *velocity;
//   double *forces;
//   double *acceleration;
//   double mass;
//   double size_time_step;
// } update_params;

// typedef struct dist_params{ 
//   int dimension; 
//   double *r1; 
//   double *r2; 
//   double *diplacement_vector;
// } dist_params;

// typedef struct random_params
// {
//   int m; 
//   int n; 
//   double a; 
//   double b; 
//   int *seed; 
//   double *r;
// }random_params;

int main ( int argc, char *argv[] )
{
  double ctime;
  pthread_barrier_init(&barrier, NULL, NUM_THREADS-1);
  pthread_mutex_init(&mutex_step, NULL);
  pthread_mutex_init(&mutex_step_print, NULL);
  pthread_mutex_init(&mutex_pe, NULL);
  pthread_mutex_init(&mutex_ke, NULL);
  pthread_mutex_init(&mutex_pe_ke, NULL);

  step_print = 0;
  step_print_index = 0;
  step_print_num = 10;
  

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
    //scanf ( "%d", &dimension );
    dimension = 2;
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
    //scanf ( "%d", &num_particles );
    num_particles = 500;
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
    //scanf ( "%d", &step_num );
    step_num = 500;
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
    //scanf ( "%lf", &size_time_step );
    size_time_step = 0.1;
  }


/*
  Allocate memory.
*/
  acceleration = ( double * ) malloc ( dimension * num_particles * sizeof ( double ) );
  forces = ( double * ) malloc ( dimension * num_particles * sizeof ( double ) );
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

  

//Pour chaque time step: 
	//printf("\n for step %d",step_num);
  pthread_t threads[NUM_THREADS];

  // long* i;
  // params* params = malloc(sizeof(params));
  // *params = (update_params) {ii, num_particles, dimension, position, velocity, acceleration};
  // params->num_particles = num_particles;

  initialize ();
  compute(-1);
  e0 = potential + kinetic;
  //printf("\n e0 = %f",e0);
  printf ( "  %8d  %14f  %14f  %14e\n", 0, potential, kinetic,
       ( potential + kinetic - e0 ) / e0 );
  step_print_index = step_print_index + 1;
  step_print = ( step_print_index * step_num ) / step_print_num;

  ctime = cpu_time ( );
  for(long i=0;  i < NUM_THREADS; i++ ){
    //printf("\n creating thread %ld",i);
    pthread_create (&threads[i], NULL, parallel, (void *)i);
  }
  for(int i=0;  i < NUM_THREADS; i++ ){
    pthread_join(threads[i], NULL);
  }
  
//Le temps d'éxecusion total:
  //printf("\n temps debut : %f", ctime);
  //printf("\n temps debut : %f", cpu_time ( ));
  ctime = cpu_time ( ) - ctime;
  printf ( "\n" );
  printf ( "  Le temps d'execution: %f seconds.\n", ctime );
/*
  Free memory.
*/
 
  pthread_mutex_destroy(&mutex_step);
  pthread_mutex_destroy(&mutex_step_print);
  pthread_mutex_destroy(&mutex_pe);
  pthread_mutex_destroy(&mutex_ke);
  pthread_mutex_destroy(&mutex_pe_ke);
  free ( acceleration );
  free ( forces );
  free ( position );
  free ( velocity );
  return 0;
}


void* parallel(void* t){
  
  
  
  long tid = (long) t;
  //printf("\n parallel function tid = %ld",tid);
  //printf("\n step_num= %d",step_num);
  for ( step = 1; step <= step_num; step++ )
  {
    //printf("\n step= %d thread= %ld",step,tid);
    pthread_mutex_lock(&mutex_step);
    int my_step = step;
    pthread_mutex_unlock(&mutex_step);
    //printf("\n chi chi tid = %ld",tid);
  	//printf("\n for step = %d", step);
      //si step == 0 on doit initialiser la position, la velocité (la vitesse) et l'acceleration de chaque particule --> check la fonction initialize en bas

    // if (( my_step == 0 )&&( tid == 0))
    // {
    //   initialize ();
    // }
    //   // sinon: on update les patricules avec les nouvelles positions, vitesses et accelerations --> check la fonction update en bas
    // else
    // {
      update(tid);
    // }
    //pthread_barrier_wait(&barrier);
    compute(tid);

  //     ke = ke * 0.5 * mass;
  // //printf("\n pe before %f",pe);
  // *(params->potential_energy) = pe;
  // //printf("\n pe after %f",pe);
  // *(params->kinetic_energy) = ke;

    


    //pthread_barrier_wait(&barrier);

    //printf("\n my_step= %d    step_print = %d",my_step,step_print);
    if ( my_step == step_print)
    {
      //printf("\n printing from thread : %ld  step = %8d step_prit = %d\n",tid, step, step_print);
      printf ( "  %8d  %14f  %14f  %14e\n", my_step, potential, kinetic,
       ( potential + kinetic - e0 ) / e0 );
      pthread_mutex_lock(&mutex_step_print);
      step_print_index = step_print_index + 1;
      step_print = ( step_print_index * step_num ) / step_print_num;
      pthread_mutex_unlock(&mutex_step_print);
      //idée : step_print = step_print + step_print_num;
    }

  }
 pthread_barrier_destroy(&barrier);
 pthread_exit(NULL);
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

void initialize ()
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
    pour dimension =2 et pour 2 paticules A et B
    le tableau est de la forme: velocityAdim1 | velocityAdim2 | velocityBdim1 | velocityBdim2
  */
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

void *update (int ii)
{
  //printf("\n update function");
  int i;
  int j;
  double rmass;
    
  rmass = 1.0 / mass;
  long start = ii * num_particles / NUM_THREADS;
  //printf("\n i = %ld \t num_part = %d \t produit = %ld \t num_thread = %d", *(params->i), num_particles, *(params->i) * num_particles, NUM_THREADS);
  long end = (ii+1) * num_particles / NUM_THREADS;
  //printf("\n Hello from update %ld start: %ld end:%ld",*(params->i),start,end);

  for ( j = start; j < end; j++ )
  {
    for ( i = 0; i < dimension; i++ )
    {
      position[i+j*dimension] = position[i+j*dimension] + velocity[i+j*dimension] * size_time_step + 0.5 * acceleration[i+j*dimension] * size_time_step * size_time_step;
      velocity[i+j*dimension] = velocity[i+j*dimension] + 0.5 * size_time_step * ( forces[i+j*dimension] * rmass + acceleration[i+j*dimension] );
      acceleration[i+j*dimension] = forces[i+j*dimension] * rmass;
    }
  }

  //return;
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


void *compute (int tid)
{
  //printf("\n compute function ");
  double d; //pour stocker la distance entre deux particules dans la boucle interne
  double d2;
  int i;
  int j;
  int k;
  double ke;
  double pe;
  //printf("\n hello from compute 3");
  double PI2 = 3.141592653589793 / 2.0;
  double rij[3];
  //printf("\n hello from compute 4");

  pe = 0.0;
  ke = 0.0;

  //printf("\n hello from compute %d", num_particles);
  //printf("\n hello from compute 5");

  long start;
    long end;
  if(tid == -1){
    start= 0;
    end= num_particles;
  }else{
    start = tid * num_particles / NUM_THREADS;
    end = (tid+1) * num_particles / NUM_THREADS;
  }
  
  //printf("\n hello from compute %ld start = %ld end = %ld", *(params->i), start, end);

//pour chaque particule
  //for ( k = 0; k < num_particles; k++ )
  for(k = start; k < end; k++ )
  {
    //printf("\n compute k = %d",k);
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
        d = dist (position+k*dimension, position+j*dimension, rij);
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
        
        //printf("\n pe before %f",pe);
        pthread_mutex_lock(&mutex_pe);
        pe = pe + 0.5 * pow ( sin ( d2 ), 2 );
        pthread_mutex_unlock(&mutex_pe);
        //printf("\n pe after %f",pe);

        for ( i = 0; i < dimension; i++ )//la force d'un particule dépend de sa distance de tous les autres  particules
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
      pthread_mutex_lock(&mutex_ke);
      ke = ke + velocity[i+k*dimension] * velocity[i+k*dimension];
      pthread_mutex_unlock(&mutex_ke);
    }
    //printf("\nk= %d pe inside %f",k,pe);
  }
  //printf("\n pe beeeeeeeefore %f",pe);

  //double tab[2] = {ke,pe};
   
  if(tid == -1 || tid==0){
    pthread_mutex_lock(&mutex_ke);
  ke = ke * 0.5 * mass;
  pthread_mutex_unlock(&mutex_ke);
    pthread_mutex_lock(&mutex_pe_ke);
    
    //printf("\n pe before %f",pe);
    
    potential = pe;
    //printf("\n pe after %f",potential);
    kinetic = ke;
    pthread_mutex_unlock(&mutex_pe_ke);
  }
  //return;
}

//Le role de cette fonction est de calculer le temps d'éxecusion d'un programme
double cpu_time ( )
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;
  //printf("\n clock = %ld",clock());

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

double dist (double r1[], double r2[], double diplacement_vector[])
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
