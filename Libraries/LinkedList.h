#include <iostream>

class Node2D{
  private:
    float count;
    //Node2D  *next;
  public:
    int x; //Not good practice but easier to implement
    int y;
    Node2D  *next;
    Node2D(int x, int y);
    Node2D();
    int addNewNode(Node2D *newNode);
    void Normalize(int totalCount);
    float getProb();
    Node2D* getNext();
    friend bool operator< (const Node2D &node1, const Node2D &node2);
    friend bool operator> (const Node2D &node1, const Node2D &node2);
    friend bool operator== (const Node2D &node1, const Node2D &node2);
    void operator= (const Node2D &node );
    static int addNewNode(Node2D **head, Node2D *newNode);
    static float getProb(Node2D *head, int posx, int posy);
};

class Node1D{   //Use for saving the marginal probabilities when needed
  private:
    float prob;
    //Node1D  *next;
  public:
    Node1D  *next;
    int z; //Not good practice but easier to implement
    Node1D(int posz, float initProb);
    Node1D();
    void addNewNode(Node1D *newNode);
    float getProb();
    Node1D* getNext();
    friend bool operator< (const Node1D &node1, const Node1D &node2);
    friend bool operator> (const Node1D &node1, const Node1D &node2);
    friend bool operator== (const Node1D &node1, const Node1D &node2);
    static void addNewNode(Node1D **head, Node1D *newNode);
};