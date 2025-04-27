#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define numberNeurons 27 //Número máximo de neuronas en el sustrato
#define numberConections 200 //Número máximo de conexiones permitidos por red
#define individuals 100 //Total de redes
#define generations 100 //Generaciones
#define mutationRate 2 //Porcentaje de mutación
#define sustitutionRate 20 //Porcentaje de sustitucion
#define elitism 1 //Elitismo posible
#define sensors 3 //Total de sensores
#define numberMaxOfInnovation 1000 //Máximo número de conexiones nuevas
#define numberMaxInitialConnections 20  //Número máximo de conexiones en las redes de la primera generación
#define initialInput 4 //Neuronas en capa de entrada máximas en la primera gen
#define initialHidden 3 //Neuronas en capa oculta máximas en la primera gen
#define inputMax 9 //Máximo de neuronas en capa de entrada
#define hiddenMax 16 //Maximo de neuronas en hidden layer
#define outputMax 2 //Máximo de salidas
#define trials 3
#define steps 6000
#define tFPossiple 4


//Variables RNA
int matrixANN[numberNeurons+sensors][numberNeurons+sensors];//Matriz que representa a la RNA y sus conexiones entre el total de nodos posibles
double outputNeurons[3][numberNeurons];//Guarda Net, FT y salida
double signalSensors[sensors];//Guarda el estado de sensores
//REVISAR Aquí, parece que están desfasadas
float Net[numberNeurons+sensors], FT[numberNeurons+sensors]; //Net y FT de la red
float biases[numberNeurons];
//Variables AGG
int parents[individuals][4];
int nodes[individuals][numberNeurons][4]; //Guarda número de nodo, bias,capa y función de transferencia de cada nodo posible
int nodesSons[individuals][numberNeurons][4]; //Guarda número de nodo, bias,capa y función de transferencia de cada nodo hijo o mutado
double fitness[generations][individuals]; //Finmtess del individuo
int listInnovationIndex[numberMaxOfInnovation][3]; //Lista para conexiones entre nodos que guarda el índice, nodo origen y destino
int theBest[elitism][numberNeurons+numberConections];//Guarda a los individuos de elitismo
int conections[individuals][numberConections][5];//Guarda nodo origen y destino, peso sináptico, habilitación e indexado
int conectionsSons[individuals][numberConections][5];//Guarda nodo origen y destino, peso sináptico, habilitación e indexado de hijos o mutaciones
int transferFunctions[numberNeurons+sensors];
int globalCounterConnections; //Cuenta la cantidad de conexiones existentes en el mapa genético. Sirve para mutación.

float fitnessIndividual;



void CleanVariables(void){
  int i, j, k;
  for(i=0;i<numberMaxOfInnovation;i++){
    listInnovationIndex[i][0]=-1;
    listInnovationIndex[i][1]=-1;
    listInnovationIndex[i][2]=-1;
  }

  for(i=0;i<individuals;i++){
    fitness[i][0]=0.0;
    if(i<numberNeurons){
      outputNeurons[0][i]=0;
      outputNeurons[1][i]=0;
      outputNeurons[2][i]=0;
    }
    if(i<sensors)
      signalSensors[i]=0.0;
    for(j=0;j<numberConections;j++){
      if(i<numberNeurons){
        if(j<numberNeurons){
          matrixANN[i][j]=0;
        }
      }
      for(k=0;k<5;k++){
        conections[i][j][k]=-1;
        conectionsSons[i][j][k]=-1;
        if(j<numberNeurons){
          if(k<4){
            nodes[i][j][k]=0;
            nodesSons[i][j][k]=0;
          }
        }
      }
    }
  }
  for(i=0;i<(numberNeurons+numberConections);i++){
    for(j=0;j<elitism;j++)
      theBest[j][i]=0;
  }
}

//Busca una conexión existente, si no está registrata, la indexa
int updateListIndex(int nodeOrigin, int nodeFinal){
  int i, index=-1;
  //Recorre el arreglo del indice de innovación
  for(i=0;i<numberMaxOfInnovation;i++){
    if(listInnovationIndex[i][0]==-1){
      listInnovationIndex[i][0]=nodeOrigin;//Inscribe nodo origen
      listInnovationIndex[i][1]=nodeFinal; //Inscribe nodo final
      listInnovationIndex[i][2]=i; //Inscribe índice de innovación
      index=i; //Guarda el índice para ser devuelto
      i=numberMaxOfInnovation+10; //Finaliza la búsqueda
    }
    if(listInnovationIndex[i][0]==nodeOrigin){
      if(listInnovationIndex[i][1]==nodeFinal){
        index=listInnovationIndex[i][2]; //Guarda el indice registrado
        i=numberMaxOfInnovation+10; //Finaliza la búsqueda
      }
    }

  }
  return index;
}
//Inicializa la población
void initPopulation(void){
  int i,j,k,l,counterConnections, input, hidden;
  FILE *archivo;
  archivo=fopen("inicia.txt", "w+");
  //Inicia el sistema
  srand(time(NULL)); //Generamos numero aleatorio en base al tiempo
  for(i=0;i<individuals;i++){
    for(j=0;j<numberNeurons;j++){
      nodes[i][j][0]=j+3;
      nodes[i][j][1]=rand() % (1000001);
      if((rand() % 10)>4)nodes[i][j][1]=nodes[i][j][1]*-1;
      if(i<9)nodes[i][j][2]=0;
      if((i>8)&&(i<17))nodes[i][j][2]=1;
      if(i>17)nodes[i][j][2]=2;
      //nodes[i][j][3]=3;
      nodes[i][j][3]=rand()%tFPossiple;
      //printf("Nodes %d %d: %d\n",i,j,nodes[i][j][3] );
      fprintf(archivo,"Individuo %d\n Nodos: %d %d %d %d\n",i,nodes[i][j][0],nodes[i][j][1],nodes[i][j][2],nodes[i][j][3]);
    }

    input=1+(rand()%(initialInput+1));
    hidden=1+(rand()%(initialHidden+1));
    //printf("input %d hidden %d\n", input, hidden);
    counterConnections=0;
    //Sensores a capa de entrada
    for(k=0;k<input;k++){
      for(l=0;l<sensors;l++){
        conections[i][counterConnections][0]=l;
        conections[i][counterConnections][1]=k+sensors;
        conections[i][counterConnections][2]=rand() % (1000001);
        if((rand() % 10)>4)conections[i][counterConnections][2]=conections[i][counterConnections][2]*-1;
        conections[i][counterConnections][3]=1;
        conections[i][counterConnections][4]=updateListIndex(conections[i][counterConnections][0], conections[i][counterConnections][1]);
        //if(i==0)printf("Connections %d %d %d: %d %d %d %d %d\n",i,k,l,conections[i][counterConnections][0],conections[i][counterConnections][1],conections[i][counterConnections][2],conections[i][counterConnections][3],conections[i][counterConnections][4]);
        fprintf(archivo,"Conexiones: %d %d %d %d %d\n",conections[i][counterConnections][0],conections[i][counterConnections][1],conections[i][counterConnections][2],conections[i][counterConnections][3],conections[i][counterConnections][4]);
        counterConnections++;

      }
    }
    //Capa de entrada a hidden
    for(k=0;k<hidden;k++){
      for(l=0;l<input;l++){
        conections[i][counterConnections][0]=l+sensors;
        conections[i][counterConnections][1]=k+inputMax+sensors;
        conections[i][counterConnections][2]=rand() % (1000001);
        if((rand() % 10)>4)conections[i][counterConnections][2]=conections[i][counterConnections][2]*-1;
        conections[i][counterConnections][3]=1;
        conections[i][counterConnections][4]=updateListIndex(conections[i][counterConnections][0], conections[i][counterConnections][1]);
        //if(i==0)printf("Connections %d %d %d: %d %d %d %d %d\n",i,k,l,conections[i][counterConnections][0],conections[i][counterConnections][1],conections[i][counterConnections][2],conections[i][counterConnections][3],conections[i][counterConnections][4]);
        fprintf(archivo,"Conexiones: %d %d %d %d %d\n",conections[i][counterConnections][0],conections[i][counterConnections][1],conections[i][counterConnections][2],conections[i][counterConnections][3],conections[i][counterConnections][4]);
        counterConnections++;
      }
    }
    //Hidden a output
    for(k=0;k<outputMax;k++){
      for(l=0;l<hidden;l++){
        conections[i][counterConnections][0]=l+inputMax+sensors;
        conections[i][counterConnections][1]=k+hiddenMax+inputMax+sensors;
        conections[i][counterConnections][2]=rand() % (1000001);
        if((rand() % 10)>4)conections[i][counterConnections][2]=conections[i][counterConnections][2]*-1;
        conections[i][counterConnections][3]=1;
        conections[i][counterConnections][4]=updateListIndex(conections[i][counterConnections][0], conections[i][counterConnections][1]);
        //if(i==0)printf("Connections %d %d %d: %d %d %d %d %d\n",i,k,l,conections[i][counterConnections][0],conections[i][counterConnections][1],conections[i][counterConnections][2],conections[i][counterConnections][3],conections[i][counterConnections][4]);
        fprintf(archivo,"Conexiones: %d %d %d %d %d\n",conections[i][counterConnections][0],conections[i][counterConnections][1],conections[i][counterConnections][2],conections[i][counterConnections][3],conections[i][counterConnections][4]);
        counterConnections++;
      }
    }
    globalCounterConnections+=counterConnections;
  }
  fclose(archivo);
}

void  setMatrix(int individualToTest){
  int i, j;
  for(i=0;i<numberNeurons+sensors;i++){
    for(j=0;j<numberNeurons+sensors;j++) {
      matrixANN[i][j]=0;
    }
    if(i<numberNeurons){
      biases[i]=nodes[individualToTest][i][1]/1000000.0;
      transferFunctions[i+sensors]=nodes[individualToTest][i][3];
    }
  }
  i=0;
  while(i<numberConections){
    if(conections[individualToTest][i][0]<0){
      i=1000;
    }
    else{
      matrixANN[conections[individualToTest][i][0]][conections[individualToTest][i][1]]=conections[individualToTest][i][2];
    }
    i++;
  }
  /*
  for(i=0;i<numberNeurons+sensors;i++){
    for(j=0;j<numberNeurons+sensors;j++) {
      if(indi==0)printf("Matrix %d %d %d: %d\n",indi, i,j, matrixANN[i][j]);
    }
  }*/
}


void upDateFT(void){
  int i;
  for(i=sensors;i<numberNeurons+sensors;i++){
    if(transferFunctions[i]==0){
      FT[i]=Net[i];
    }
    if(transferFunctions[i]==1)
    {
      FT[i]=(exp(Net[i])-exp(-1*Net[i]))/(exp(Net[i])+exp(-1*Net[i]));
      //if(FT[i]<0)FT[i]=-1;
      //else FT[i]=1;
    }
    if(transferFunctions[i]==2){
      FT[i]=1/1+exp(-1*Net[i]);
      //if(FT[i]<0.5)FT[i]=0;
      //else FT[i]=1;
    }
    if(transferFunctions[i]==3){
      FT[i]=Net[i];
      if(FT[i]<-1.0)FT[i]=-1.0;
      if(FT[i]>1.0)FT[i]=1.0;
      //FT[i]=1/1+exp(-1*Net[i]);
      //if(FT[i]<0.5)FT[i]=-1;
      //if(FT[i]==0.5)FT[i]=0;
      //if(FT[i]>0.5)FT[i]=1;
    }
    if(transferFunctions[i]==4)
    {
      FT[i]=(exp(Net[i])-exp(-1*Net[i]))/(exp(Net[i])+exp(-1*Net[i]));
      if(FT[i]<0)FT[i]=-1;
      if(FT[i]==0)FT[i]=0;
      if(FT[i]>0)FT[i]=1;
    }
    if(transferFunctions[i]==5)
    {
      FT[i]=(exp(Net[i])-exp(-1*Net[i]))/(exp(Net[i])+exp(-1*Net[i]));
    }
    if(transferFunctions[i]==6){
      FT[i]=1/1+exp(-1*Net[i]);
    }

  }
}

//Actualiza motores
void motorsUpDate(void){
  int i=sensors+inputMax+hiddenMax;
  leftSpeed=FT[i];
  i++;
  rightSpeed=FT[i];
  //if(rightSpeed>0||leftSpeed>0)printf("Salida: %f %f\n",leftSpeed,rightSpeed);
  int obstacles=0;
  //printf("R:%f L:%f %d\n",leftSpeed,rightSpeed, obstacles);
  //Contabiliza fitness
  for(i=0;i<sensors;i++){
    //fitnessIndividual+=(1000-signalSensors[i]);
    if(signalSensors[i]<990)obstacles++;
    //printf("Sensor %d: %f",i, signalSensors[i]);
    //fitnessIndividual+=signalSensors[i];
  }
  //printf("R:%f L:%f %d\n",leftSpeed,rightSpeed, obstacles);

  //if((rightSpeed==0)&&(leftSpeed==0))fitnessIndividual=fitnessIndividual-1;
  //else fitnessIndividual=fitnessIndividual+(1-sqrt(pow((rightSpeed-leftSpeed),2)));
  //fitnessIndividual=fitnessIndividual+sqrt(pow((rightSpeed+leftSpeed),2));
  if((rightSpeed!=0)&&(leftSpeed!=0)){
    //printf("Fitness vel dif 0: %f\n",fitnessIndividual);
    fitnessIndividual=fitnessIndividual+(3-obstacles);
    fitnessIndividual=fitnessIndividual+(1-sqrt(pow((rightSpeed-leftSpeed),2)));
    //printf("Fitness vel dif 0: %f, Obstacles %d, %f %f\n",fitnessIndividual,obstacles,rightSpeed,leftSpeed);
  }
  else{
    fitnessIndividual=fitnessIndividual-5;
    //printf("Fitness vel= 0: %f\n",fitnessIndividual);
  }

  //printf("Velocidades = %f %f %d %f %f %f\n", leftSpeed, rightSpeed, obstacles, signalSensors[0],signalSensors[1],signalSensors[2]);
}

//Envia el mensaje de finalización al supervisor
void Send(void){
   char message[128]={"1"};
   wb_emitter_send(rfe, message, strlen(message) + 1);
}

//Computa el estado de la RNA
void ANNCompute(void){
  int i,j; //Índices
  float inputs[numberNeurons+sensors]; //Entradas de la red
  for(i=0;i<(numberNeurons+sensors);i++){
    inputs[i]=FT[i]; //Limpia entradas
    Net[i]=0.0f; //NET
  }
  for(i=0;i<sensors;i++){
    //inputs[i]=1.0-(signalSensors[i]/1000.0); //Guarda sensores
    inputs[i]=signalSensors[i]/1000.0;
  }
  //Calcula Input
  for(i=sensors;i<(sensors+inputMax);i++){
    for(j=0;j<(numberNeurons+sensors);j++){
      Net[i]+=(matrixANN[j][i]/1000000.0*inputs[j]);
      //printf("Net %d: %f, Peso: %d, Input: %f\n", i, Net[i],matrixANN[j][i],inputs[j]);
    }
    Net[i]+=biases[i-sensors];
    //printf("Net: %f, Bias: %f\nFin neurona\n", Net[i],biases[i-sensors]);
  }
  upDateFT();
  for(i=sensors;i<sensors+inputMax;i++){
    inputs[i]=FT[i]; //las salidas calculadas
  }
  //Hidden
  for(i=(sensors+inputMax);i<(sensors+inputMax+hiddenMax);i++){
    for(j=0;j<(numberNeurons+sensors);j++){
      Net[i]+=(matrixANN[j][i]/1000000.0)*inputs[j];
      //printf("Net: %f, Peso: %d, Input: %f\n", Net[i],matrixANN[j][i],inputs[j]);
    }
    Net[i]+=biases[i-sensors];
    //printf("Net: %f, Bias: %f\nFin neurona\n", Net[i],biases[i-sensors]);
  }
  upDateFT();
  for(i=sensors+inputMax;i<sensors+inputMax+hiddenMax;i++){
    inputs[i]=FT[i]; //las salidas calculadas
  }
  //Output
  for(i=(sensors+inputMax+hiddenMax);i<(sensors+inputMax+hiddenMax+outputMax);i++){
    for(j=0;j<(numberNeurons+sensors);j++){
      Net[i]+=(matrixANN[j][i]/1000000.0)*inputs[j];
      //printf("Net: %f, Peso: %d, Input: %f\n", Net[i],matrixANN[j][i],inputs[j]);
    }
    Net[i]+=biases[i-sensors];
    //printf("Net: %f, Bias: %f\nFin neurona\n", Net[i],biases[i-sensors]);
  }
  upDateFT();
  //printf("Salidas\n");
  for(i=0;i<numberNeurons+sensors;i++){
    inputs[i]=FT[i]; //las salidas calculadas
    //printf("%d %f %f %f %d\n",i,inputs[i], Net[i],FT[i],transferFunctions[i]);
  }
  motorsUpDate();
}

void Test(void){
  int i,j,k;
  fitnessIndividual=0.0f;
  for(k=0;k<(numberNeurons+sensors);k++){
      FT[k]=0.0f; //FT
  }
  for(i=0;i<trials;i++){
    for(j=0;j<steps;j++){
      SenseStep();
      ANNCompute();
    }
  printf("Fin Trial\n");
  }


}

//Selecciona los padres potenciales
void Selection(int generationIndex){
  int i, a, b;
  //FILE *archivo;
  //archivo=fopen("seleccion.txt", "w+");
  for(i=0;i<individuals;i++){
     //parent A
     a=rand()%(individuals);
     b=rand()%(individuals);
     if(b==a){
       if(b<99)b++;
       else b--;
     }
     //fprintf(archivo,"A1: %d, %f B1: %d, %f\n", a, fitness[generationIndex][a],  b,fitness[generationIndex][b]);
     if(fitness[generationIndex][a]>fitness[generationIndex][b])parents[i][0]=a;
     else parents[i][0]=b;
     //Parent B
     a=rand()%(individuals);
     b=rand()%(individuals);
     if(b==a){
       if(b<99)b++;
       else b--;
     }
     //fprintf(archivo,"A2: %d, %f B2: %d, %f\n", a,  fitness[generationIndex][a], b,fitness[generationIndex][b]);
     if(fitness[generationIndex][a]>fitness[generationIndex][b])parents[i][1]=a;
     else parents[i][1]=b;
  }
  //fclose(archivo);
}

//Realiza la cruza de los padres seleccionados
void CrossOver(int generationIndex){
  //FILE *archivo;
  //archivo=fopen("crossover.txt", "w+");
  printf("Entra a cruza\n");
  int A, B,i,j,k,crossPoint, flag; //La bandera sirve para indicar si se enconraron dos conexiones iguales
  for(i=0;i<individuals;i++){
    //El padre A es aquel con el fitness más alto
    A=parents[i][0];
    B=parents[i][1];
    if(fitness[generationIndex][B]>fitness[generationIndex][A]){
      A=parents[i][1];
      B=parents[i][0];
    }
    //Compara si el índice de innovación es igual entre un par de nodos
    //Si esto ocurre, copia el gen de uno de los padres al hijo.
    //Si no aparece, se copia del padre con mayor fitness
    for(j=0;j<numberConections;j++){
      flag=0;
      for(k=0;k<numberConections;k++){
        if(conections[A][j][4]==conections[B][k][4]){
          crossPoint=rand() % (10);
          flag=1;//Nodo encontrado en ambos padres
          //Padre B
          if(crossPoint>5){
            conectionsSons[i][j][0]=conections[B][k][0];
            conectionsSons[i][j][1]=conections[B][k][1];
            conectionsSons[i][j][2]=conections[B][k][2];
            conectionsSons[i][j][3]=conections[B][k][3];
            conectionsSons[i][j][4]=conections[B][k][4];
            //fprintf(archivo,"B: %d, %d, %d, %d, %d\n", conections[B][k][0], conections[B][k][1], conections[B][k][2], conections[B][k][3],conections[B][k][4]);
          }
          //Padre A
          else{
            conectionsSons[i][j][0]=conections[A][j][0];
            conectionsSons[i][j][1]=conections[A][j][1];
            conectionsSons[i][j][2]=conections[A][j][2];
            conectionsSons[i][j][3]=conections[A][j][3];
            conectionsSons[i][j][4]=conections[A][j][4];
            //fprintf(archivo,"A: %d, %d, %d, %d, %d\n", conections[A][j][0], conections[A][j][1], conections[A][j][2], conections[A][j][3],conections[A][j][4]);
          }
          k+=numberNeurons; //Sale del ciclo
        }
      }
      //Verifica si encontró conexiones iguales. Si no las encontró, copia el gen del padre más apto
      if(flag==0){
        conectionsSons[i][j][0]=conections[A][j][0];
        conectionsSons[i][j][1]=conections[A][j][1];
        conectionsSons[i][j][2]=conections[A][j][2];
        conectionsSons[i][j][3]=conections[A][j][3];
        conectionsSons[i][j][4]=conections[A][j][4];
        //fprintf(archivo,"A diferente: %d, %d, %d, %d, %d\n", conections[A][j][0], conections[A][j][1], conections[A][j][2], conections[A][j][3],conections[A][j][4]);
      }
    }
    //Ahora recorre el arreglo de nodos sorteando el gen padre
    for(j=0;j<numberNeurons;j++){
      crossPoint=rand() % (10);
     //Padre B
      if(crossPoint>5){
        nodesSons[i][j][0]=nodes[B][j][0];
        nodesSons[i][j][1]=nodes[B][j][1];
        nodesSons[i][j][2]=nodes[B][j][2];
        nodesSons[i][j][3]=nodes[B][j][3];
        //fprintf(archivo,"B: %d, %d, %d, %d\n", nodes[B][j][0], nodes[B][j][1], nodes[B][j][2], nodes[B][j][3]);
      }
      //Padre B
      else{
        nodesSons[i][j][0]=nodes[A][j][0];
        nodesSons[i][j][1]=nodes[A][j][1];
        nodesSons[i][j][2]=nodes[A][j][2];
        nodesSons[i][j][3]=nodes[A][j][3];
        //fprintf(archivo,"A: %d, %d, %d, %d\n", nodes[A][j][0], nodes[A][j][1], nodes[A][j][2], nodes[A][j][3]);
      }
    }
  }
  //fclose(archivo);
  printf("Termina cruza\n");

}

//Selecciona los individuos de la nueva generación
void Sustitution(void){
  float limit=individuals*(sustitutionRate*.01); //Establece el límite de padres que se pueden conservar
  int counter=0; //contador de padres conservados;
  int selection,i,j; //Para el aleatorio de sustitución
  //FILE *archivo;
  //archivo=fopen("substitution.txt", "w+");
  for(i=0;i<individuals;i++){
    if(counter<=limit){
      selection=rand() % (10);
      //fprintf(archivo,"Ind %d Sel: %d, Count: %d\n",i, selection, counter);
      if(selection>5){
        counter++; //Incrementa el contador de sustitución
        //Copia al padre intacto y elimina al hijo
        for(j=0;j<numberNeurons;j++){
          nodesSons[i][j][0]=nodes[i][j][0];
          nodesSons[i][j][1]=nodes[i][j][1];
          nodesSons[i][j][2]=nodes[i][j][2];
          nodesSons[i][j][3]=nodes[i][j][3];
        }
        for(j=0;j<numberConections;j++){
          conectionsSons[i][j][0]=conections[i][j][0];
          conectionsSons[i][j][1]=conections[i][j][1];
          conectionsSons[i][j][2]=conections[i][j][2];
          conectionsSons[i][j][3]=conections[i][j][3];
          conectionsSons[i][j][4]=conections[i][j][4];
        }
      }
    }
  }
  //fclose(archivo);
}

void Mutation(void){
  int i,j,solution,origin,end, flag, onoff;
  //De conexiones
  float mutationLimit=globalCounterConnections*(mutationRate/100.0);
  //printf("MutationLimit %f globalCounterConnections %d \n",mutationLimit, globalCounterConnections);
  //FILE *archivo;
  //archivo=fopen("mutation.txt", "w+");
  for(i=0;i<mutationLimit;i++){
    solution=rand() % (individuals);
    onoff = rand() % 10;
    //Agrega conexiones
    if(onoff<5){
      origin=rand() % (numberNeurons);
      end=rand() % (numberNeurons);
      flag=0;
      j=0;
      //fprintf(archivo,"Sorteo: ind %d inicio %d fin %d", solution, origin, end);
      //Busca si la conexión sorteada se encuentra en la red
      while(flag==0){
        if(conectionsSons[solution][j][0]>-1){
          if(conectionsSons[solution][j][0]==origin){
            if(conectionsSons[solution][j][1]==end){
              flag=1; //Si la encuentra, termina la búsqueda sin hacer nada
              //fprintf(archivo, "Mutación incluida en j %d\n", j);
            }
          }
        }
        else{ //Si no la encuentra y llega al final de la lista
          conectionsSons[solution][j][4]=updateListIndex(origin, end);//Inscribe la conexión mutada
          conectionsSons[solution][j][0]=origin;
          conectionsSons[solution][j][1]=end;
          conectionsSons[solution][j][3]=1;
          conectionsSons[solution][j][2]=rand() % (1000001);
          if((rand() % 10)>4)conectionsSons[solution][j][2]=conectionsSons[solution][j][2]*-1;
          //fprintf(archivo,"Conection mutada a individuo %d: %d %d %d %d %d\n", solution, conectionsSons[solution][j][0],conectionsSons[solution][j][1],conectionsSons[solution][j][2],conectionsSons[solution][j][3],conectionsSons[solution][j][4]);
          flag=1; //Finaliza el cilco
          globalCounterConnections=globalCounterConnections+1;
        }
        j++; //Incrementa el contador del buscador
      }
    }
    //Retira conexiones
    else{
      j = rand() % numberConections;
      conectionsSons[solution][j][0]=-1;
      conectionsSons[solution][j][1]=-1;
      conectionsSons[solution][j][2]=-1;
      conectionsSons[solution][j][3]=-1;
      conectionsSons[solution][j][4]=-1;
      globalCounterConnections=globalCounterConnections-1;
    }
  }
  //De FT
  mutationLimit=numberNeurons*individuals*(mutationRate/100.0);
  int tF,neuron; //Función de transferencia y neurona
  for(i=0;i<mutationLimit;i++){
    solution=rand() % (100);
    tF=rand()%tFPossiple;
    neuron=rand() % (numberNeurons);
    nodesSons[solution][neuron][3]=tF;
    //fprintf(archivo,"Sorteo TF: ind %d tf %d neuron %d\n", solution, tF, neuron);
  }
  //fclose(archivo);
}

void Elitism(int gene){
  int i,mayor;
  mayor=0;
  //Busca el mayor fitness
  for(i=1;i<individuals;i++){
    if(fitness[gene][i]>fitness[gene][mayor])mayor=i;
  }
  //Guarda el mejor de los individuos
  //Nodos
  for(i=0;i<numberNeurons;i++){
    nodesSons[mayor][i][0]=nodes[mayor][i][0];
    nodesSons[mayor][i][1]=nodes[mayor][i][1];
    nodesSons[mayor][i][2]=nodes[mayor][i][2];
    nodesSons[mayor][i][3]=nodes[mayor][i][3];
  }
  //Conexiones
  for(i=0;i<numberConections;i++){
    conectionsSons[mayor][i][0]=conections[mayor][i][0];
    conectionsSons[mayor][i][1]=conections[mayor][i][1];
    conectionsSons[mayor][i][2]=conections[mayor][i][2];
    conectionsSons[mayor][i][3]=conections[mayor][i][3];
    conectionsSons[mayor][i][4]=conections[mayor][i][4];
  }
  //printf("Mayor generacion %d %d: %f", mayor, gene, fitness[gene][mayor]);
}

void NewGeneration(void){
  int i,j;
  for(i=0;i<individuals;i++){
    for(j=0;j<numberNeurons;j++){
      nodes[i][j][0]=nodesSons[i][j][0];
      nodes[i][j][1]=nodesSons[i][j][1];
      nodes[i][j][2]=nodesSons[i][j][2];
      nodes[i][j][3]=nodesSons[i][j][3];
      nodesSons[i][j][0]=0;
      nodesSons[i][j][1]=0;
      nodesSons[i][j][2]=0;
      nodesSons[i][j][3]=0;
    }
   for(j=0;j<numberConections;j++){
     conections[i][j][0]=conectionsSons[i][j][0];
     conections[i][j][1]=conectionsSons[i][j][1];
     conections[i][j][2]=conectionsSons[i][j][2];
     conections[i][j][3]=conectionsSons[i][j][3];
     conections[i][j][4]=conectionsSons[i][j][4];
     conectionsSons[i][j][0]=-1;
     conectionsSons[i][j][1]=-1;
     conectionsSons[i][j][2]=-1;
     conectionsSons[i][j][3]=-1;
     conectionsSons[i][j][4]=-1;
   }
  }
}

void PrintGen(int gen){
  int i,j;
  FILE *archivo, *hijos;
  archivo=fopen("resultados.txt", "w+");
  hijos=fopen("hijos.txt","w+");
  fprintf(archivo,"Generacion %d\n",gen);
  for(i=0;i<individuals;i++){
    fprintf(archivo,"Individuo: %d ", i);
    fprintf(hijos,"Individuo: %d ", i);
    for(j=0;j<numberNeurons;j++){
      fprintf(archivo," Nodo [%d]: %d, %d %d %d\n", j, nodes[i][j][0],nodes[i][j][1],nodes[i][j][2],nodes[i][j][3]);
      fprintf(hijos," Nodo [%d]: %d, %d %d %d\n", j, nodesSons[i][j][0],nodesSons[i][j][1],nodesSons[i][j][2],nodesSons[i][j][3]);
    }
    fprintf(archivo,"\n");
    fprintf(hijos,"\n");
    for(j=0;j<numberConections;j++){
      if(conections[i][j][0]!=-1){
        fprintf(archivo," Conexión [%d]: %d, %d %d %d %d\n", j, conections[i][j][0],conections[i][j][1],conections[i][j][2],conections[i][j][3],conections[i][j][4]);
      }
      else{
        j+=numberConections;
      }
    }
    for(j=0;j<numberConections;j++){
      if(conectionsSons[i][j][0]!=-1){
        fprintf(hijos," Conexión [%d]: %d, %d %d %d %d\n", j, conectionsSons[i][j][0],conectionsSons[i][j][1],conectionsSons[i][j][2],conectionsSons[i][j][3],conectionsSons[i][j][4]);
      }
      else{
        j+=numberConections;
      }
    }
    fprintf(archivo,"\nFitness [%d]: %f\n", i, fitness[gen][i]);
  }
  fclose(archivo);
  fclose(hijos);
}
//Imprime fitness
void PrintFitness(void){
  int i,j;
  FILE *archivo;
  archivo=fopen("fitness.txt", "w+");
  for(j=0;j<generations;j++){
    for(i=0;i<individuals;i++){
      fprintf(archivo," Generación [%d] Individuo [%d]: %f\n", j,i,fitness[j][i]);
    }
  }
  fclose(archivo);
}

void PrintVariables(void){
  int i,j;
  for(i=0;i<3;i++){
    printf("signalSensors %d %d: %f\n", i,j,signalSensors[i]);
    for(j=0;j<27;j++){
      printf("outputNerons %d %d: %f\n", i,j,outputNeurons[i][j]);
    }
  }

  for(i=0;i<27;i++){
    for(j=0;j<27;j++){
      if(j<4){
        printf("node %d %d: %d\n", i,j,nodes[0][i][j]);
        printf("nodeSons %d %d: %d\n", i,j,nodesSons[0][i][j]);
      }
        printf("Matrix %d %d: %d\n", i,j,matrixANN[i][j]);
    }
  }
  for(i=0;i<100;i++){
    for(j=0;j<100;j++){
      if(j<3)printf("listInnovation %d %d: %d\n", i,j,listInnovationIndex[i][j]);
      printf("Fitness %d %d: %f\n", i,j,fitness[i][j]);
    }
  }

  for(j=0;j<227;j++){
      printf("theBest %d, %d: %d\n", i,j,theBest[0][j]);
  }
}

void PrintConnections(void){
  for(int i=0;i<100;i++){
    for(int j=0;j<200;j++){
      for(int k=0; k<5;k++){
        printf("conections %d %d %d: %d\n", i,j,k, conections[i][j][k]);
        printf("conectionsSons %d %d %d: %d\n", i,j,k, conectionsSons[i][j][k]);
      }
    }
  }
}

void Final(void){
    PrintFitness();
    printf("\nFinal de la evolución");
    int max=0;
    for(int i=1;i<individuals;i++){
      if(fitness[generations-1][i]>fitness[generations-1][max])
        max=i;
    }
    setMatrix(max);
}


int main(int argc, char *argv[]) {

  int individualCounter,generationCounter;

  CleanVariables();
  //PrintVariables();
  //PrintConnections();
  globalCounterConnections=0;
  initPopulation();
  while (wb_robot_step(time_step) != -1) {
    for(generationCounter=0;generationCounter<generations;generationCounter++){
      printf("Inicia generacion %d", generationCounter);
      for(individualCounter=0;individualCounter<individuals;individualCounter++){
        setMatrix(individualCounter);
        Test();
        fitness[generationCounter][individualCounter]=fitnessIndividual;
        printf("%d\n",individualCounter);
        wb_robot_step(time_step);
      }
      //operatorsApplication(generationCounter);
      printf("Inicia seleccion\n");
      Selection(generationCounter);
      printf("Inicia cruza\n");
      //wb_robot_step(time_step);
      CrossOver(generationCounter);
      printf("Inicia sustitucion\n");
      //wb_robot_step(time_step);
      Sustitution();
      printf("Inicia mutacion\n");
      Mutation();
      Elitism(generationCounter);
      PrintGen(generationCounter);
      NewGeneration();
      wb_robot_step(time_step);

    }
    Final();
    int i=0;
    while(i==0){
      SenseStep();
      ANNCompute();
      wb_robot_step(time_step);
    }

  }
  //Sense();
  wb_robot_cleanup();

  return 0;

}
