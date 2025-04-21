/*-------------------------------------------------------
                    <copyright notice>

Created by Dongyue Bai
date: 18.01.2019
In this work, adjacency matrix representation is used for graph representation. 

---------------------------------------------------------*/
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <random>
#include <limits>


//creat a class for graph
class Graph{
public:

bool *ShortestPathTree; //visited vertices will be marked as true, unvisited vertices will be marked as false
double *dist;  // used to store the distance from source to i-th element

//constructors
Graph();

Graph(int _size,double _density, double _distance_range_low, double _distance_range_high): 
m_size(_size), m_density(_density), m_distance_range_low(_distance_range_low),m_distance_range_high(_distance_range_high){};  


// function that draws a randomly generated 2D graph *m_graph[] and put edges in the graph by a density
void f_draw_graph(); 

//function that prints the generated graph
void f_print_graph(); 

//function that checks if the graph is connected
bool f_is_connected(); 

// function that implments dijkstras algorithm
void dijkstra();


//destructor
~Graph(){
    delete[] m_graph;
};


private:
double **m_graph;
int m_size;
double m_density;
double m_distance_range_low;
double m_distance_range_high;
};


void Graph::f_draw_graph()
{
//create a graph
m_graph= new double*[m_size];
for(int i=0;i<m_size;++i)  {m_graph[i]=new double[m_size];};    

 //seed rand()
srand(time(0));

//put edges in the graph
for(int i=0;i<m_size;++i)
{
   for(int j=0;j<m_size;++j)
   {
        if(i==j) 
        {
            m_graph[i][j]=0;  
        }
        else 
        {
            
            bool edge=((rand() % 100) < m_density*100);  // the probability of edge to be 1 is the value of density
            
            double f = (double)rand()/RAND_MAX;
            double distance= m_distance_range_low + f * (m_distance_range_high - m_distance_range_low);
            m_graph[i][j]=m_graph[j][i]=(static_cast<double> (edge))*distance;
        }
    }     
}
}

void Graph::f_print_graph()
{
    for(int i=0;i<m_size;++i){
        for(int j=0;j<m_size;++j){
            std::cout<<m_graph[i][j]<<"\t";
        }
        std::cout<<std::endl;
    }
} 


//check if the graph is connected
bool Graph::f_is_connected()  
{
    int old_size, c_size=0;
    bool* close=new bool[m_size];
    bool* open=new bool[m_size];
    for(int i=0;i<m_size;++i)
    {
        close[i]=open[i]=false; //initialization, both close set and open set are empty
    }
    open[0]=true;


    while(c_size<m_size)
    {
        for(int i=0;i<m_size;++i)
        {
            old_size=c_size;
            if(open[i]&&(close[i]==false))
            {   close[i]=true; 
                c_size++;    //add to close set 
                for(int j=0;j<m_size;++j)
                {
                    open[j]=open[j]||(m_graph[i][j]>0);   //add to open set
                }
            }
        }
        if(c_size==m_size)    {return true;}  //if all nodes in close set, connected
        if(c_size== old_size) {return false;}
       
    }
}




// Function that implements Dijkstra's single source shortest path algorithm 
// for a graph represented using adjacency matrix representation 
void Graph::dijkstra()
{
    bool* ShortestPathTree=new bool[m_size];
    double* dist=new double[m_size];

     // Initialize all distances as INTMAX and ShortestPathTree[] as false 
     for (int i = 0; i <m_size; i++)
     {
        dist[i] = std::numeric_limits<int>::max();
        ShortestPathTree[i] = false; 
     }
     
     // Distance of source vertex from itself is always 0 
     dist[0] = 0; 

     // Find shortest path for all vertices 
     for (int count = 0; count < m_size-1; ++count) 
     { 

       // Pick the minimum distance vertex from the set of vertices not yet visited. 
        double min = std::numeric_limits<int>::max();
        int min_index; 
        for (int i = 0; i < m_size; ++i) 
        {
            if (!ShortestPathTree[i] && dist[i]<=min )
                {
                    min = dist[i]; 
                    min_index = i;                     
                }             
        }

       int u = min_index; 
   
       // Mark the picked vertex as visited 
       ShortestPathTree[u] = true; 
   
       // Update dist value of the adjacent vertices of the picked vertex. 
       for (int i = 0; i < m_size; ++i)
       {   
         // Update dist[i] only if is not in ShortestPathTree, there is an edge from  
         // u to i, and total weight of path from source to  i through u is  
         // smaller than current value of dist[i] 
         if (!ShortestPathTree[i] && m_graph[u][i] && dist[u]+m_graph[u][i] < dist[i])
            { 
            dist[i] = dist[u] + m_graph[u][i]; 
            }
        }
     } 

    //print out the result
    double sum, aver_dist;
    std::cout<<"Vertex Distance from Source"<<std::endl; 
    for (int i = 0; i < m_size; ++i)
    { 
        std::cout<<i<<" "<<dist[i]<<std::endl; 
        sum += dist[i];
    }
    aver_dist=sum/49;

    std::cout<<"the average distance is "<<aver_dist<<std::endl;
    
   
} 



int main(){
    std::cout<< "draw a graph:20% of 50 nodes with a distance range of 1.0 to 10.0"<<std::endl; 
    //generate graph
    Graph obj_graph(50,0.2,1,10);
    obj_graph.f_draw_graph();
    obj_graph.f_print_graph();

    //check if the graph is connected
    bool connection_status=obj_graph.f_is_connected();
    if (connection_status==true){
        std::cout<<"The generated graph is connected"<<std::endl;
        }
    else{std::cout<<"The generated graph is not connected"<<std::endl;}

    if(connection_status==true)
    {
        obj_graph.dijkstra();
    }
    
    
    std::cout<< "draw a graph: 40% of 50 nodes with a distance range of 1.0 to 10.0"<<std::endl; 
    //generate graph
    Graph obj2_graph(50,0.4,1,10);
    obj2_graph.f_draw_graph();
    obj2_graph.f_print_graph();

    //check if the graph is connected
    connection_status=obj2_graph.f_is_connected();
    if (connection_status==true){
        std::cout<<"The generated graph is connected"<<std::endl;
        }
    else{std::cout<<"The generated graph is not connected"<<std::endl;}

    if(connection_status==true)
    {
        obj2_graph.dijkstra();
    }

}