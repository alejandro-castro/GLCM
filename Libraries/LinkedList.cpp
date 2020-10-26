#include "LinkedList.h"

Node2D::Node2D (int posx, int posy): x(posx), y(posy), next(NULL), count(1.0){}
Node2D::Node2D(){}

int Node2D::addNewNode(Node2D *newNode){
  if (*this == *newNode){
    this->count = this->count + 1;
    return 0;
  }
  else if (this->next == NULL){
    this->next = newNode;
    return 1;
  }
  else{
    Node2D *nextNode = this->next;
    if (*newNode < *nextNode){
      this->next = newNode;
      newNode->next = nextNode;
      return 1;
    }
    else{
      return (nextNode->addNewNode(newNode));
    } 
  }
}

bool operator< (const Node2D &node1, const Node2D &node2){
  return (node1.x < node2.x) || ((node1.x == node2.x) && (node1.y < node2.y));
}

bool operator> (const Node2D &node1, const Node2D &node2){
  return (node1.x > node2.x) || ((node1.x == node2.x) && (node1.y > node2.y));
}

bool operator== (const Node2D &node1, const Node2D &node2){
  return (node1.x == node2.x) && (node1.y == node2.y);
}
void Node2D::operator= (const Node2D &node ) { 
this->x = node.x;
this->y = node.y;
this->next = node.next;
this->count = node.count;
}

void Node2D::Normalize(int totalCount){
  this->count = this->count / totalCount;
  Node2D *temp = this->next;
  if (temp != NULL)  temp->Normalize(totalCount);
}

float Node2D::getProb(){  return this->count;}

int Node2D::addNewNode(Node2D **head, Node2D *newNode){
    if (*head == NULL){
      *head = newNode;
      return 1;
    }
    else{
        if (*newNode < **head){
          Node2D *temp = *head;
          *head = newNode; // The new node becomes the new head
          (*head) -> addNewNode(temp);
          return 1;
        }
        else{
          return ((*head) -> addNewNode(newNode));
        }
    } 

}

float Node2D::getProb(Node2D *head, int posx, int posy){
    Node2D *temp = head;
    Node2D newNode = Node2D(posx, posy);
    float prob = 0.0;

    while((temp != NULL)&&!(newNode<*temp)){   
        if(newNode == *temp){
            prob = temp->getProb();
            break;
        }
        else{
            temp = temp->next;
        }   
    }

    return prob;
}

Node2D* Node2D::getNext(){ return this->next;}






Node1D::Node1D (int posz, float initProb): z(posz), next(NULL), prob(initProb){}
Node1D::Node1D(){}

void Node1D::addNewNode(Node1D *newNode){
  if (*this == *newNode){
    this->prob = this->prob + newNode->getProb();
  }
  else if (this->next == NULL){
    this->next = newNode;
  }
  else{
    Node1D *nextNode = this->next;
    if (*newNode < *nextNode){
      this->next = newNode;
      newNode->next = nextNode;
    }
    else{
      nextNode->addNewNode(newNode);
    } 
  }
}

bool operator< (const Node1D &node1, const Node1D &node2){
  return node1.z < node2.z;
}

bool operator> (const Node1D &node1, const Node1D &node2){
  return node1.z > node2.z;
}

bool operator== (const Node1D &node1, const Node1D &node2){
  return node1.z == node2.z;
}

float Node1D::getProb(){  return this->prob;}


void Node1D::addNewNode(Node1D **head, Node1D *newNode){
    if (*head == NULL){
      *head = newNode;
    }
    else{
        if (*newNode < **head){
          Node1D *temp = *head;
          *head = newNode; // The new node becomes the new head
          (*head) -> addNewNode(temp);
        }
        else{
          (*head) -> addNewNode(newNode);
        }
    } 

}

Node1D* Node1D::getNext(){ return this->next;}


