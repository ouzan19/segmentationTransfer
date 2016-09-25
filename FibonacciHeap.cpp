
// ====================================================================
//
//	A network optimization algorithm using Fibonacci Heap
//
//	Written by: Max Winkler
//
// ====================================================================

#include "FibonacciHeap.h"
#include <vector>

using namespace Dijsktra;
// =========================================================================
//	Implementation of class Edge
// =========================================================================
Edge::Edge(Node* tail, Node* head, double length)
{
	this->tail = tail;
	this->head = head;
	this->length = length;
}

// =========================================================================
//	Implementation of class Node
// =========================================================================

Node::Node(int data, double key)
{
	this->key = key;
	this->data = data;
	parent = NULL;
	children = NULL;
	leftSibling = NULL;
	rightSibling = NULL;
	pred = NULL;
	rank = 0;
	state = UNLABELED;
}

Node::Node()
{
	parent = NULL;
	children = NULL;
	leftSibling = NULL;
	rightSibling = NULL;
	pred = NULL;
	rank = 0;
	state = UNLABELED;
}

bool Node::addChild(Node *node)
{
	if(children != NULL)
		children->addSibling(node);
	else
	{
		children = node;
		node->parent = this;
		rank = 1;
	}

	return true;
}

bool Node::addSibling(Node *node)
{
	Node* temp = rightMostSibling();
	if(temp == NULL)
		return false;

	temp->rightSibling = node;
	node->leftSibling = temp;
	node->parent = this->parent;
	node->rightSibling = NULL;

	if(parent)
		parent->rank++;

	return true;
}

bool Node::remove()
{
	if(parent)
	{
		parent->rank--;
		if(leftSibling)
			parent->children = leftSibling;
		else if(rightSibling)
			parent->children = rightSibling;
		else
			parent->children = NULL;
	}	
	
	if(leftSibling)
		leftSibling->rightSibling = rightSibling;
	if(rightSibling)
		rightSibling->leftSibling = leftSibling;
	
	leftSibling = NULL;
	rightSibling = NULL;
	parent = NULL;

	return true;
}

void Node::addIncomingEdge(Edge * edge)
{
	incomingEdges.push_back(edge);
}

void Node::addOutgoingEdge(Edge * edge)
{
	outgoingEdges.push_back(edge);
}


Node* Node::leftMostSibling()
{
	if(this == NULL)
		return NULL;

	Node* temp = this;
	while(temp->leftSibling)
		temp = temp->leftSibling;
	return temp;
}

Node* Node::rightMostSibling()
{
	if(this == NULL)
		return NULL;

	Node* temp = this;
	while(temp->rightSibling)
		temp = temp->rightSibling;
	return temp;
}

// =========================================================================
//	Implementation of class Fibonacci Heap
// =========================================================================


FibonacciHeap::FibonacciHeap()
{
	minRoot = NULL;
}

FibonacciHeap::FibonacciHeap(Node *root)
{
	this->minRoot = root;
	minRoot->parent = NULL;
	minRoot->children = NULL;
	minRoot->leftSibling = NULL;
	minRoot->rightSibling = NULL;
	minRoot->rank = 0;
}

FibonacciHeap::~FibonacciHeap()
{
	delete[] rootListByRank;
}

bool FibonacciHeap::isEmpty()
{
	return (minRoot == NULL);
}

bool FibonacciHeap::insertVertex(Node * node)
{
	if(node == NULL)
		return false;

	if(minRoot == NULL)
		minRoot = node;
	else
	{
		minRoot->addSibling(node);
		if(minRoot->key > node->key)
			minRoot = node;
	}
	return true;
}

Node* FibonacciHeap::findMin()
{
	return minRoot;
}

Node* FibonacciHeap::deleteMin()
{
	Node *temp = minRoot->children->leftMostSibling();
	Node *nextTemp = NULL;

	// Adding Children to root list		
	while(temp != NULL)
	{
		nextTemp = temp->rightSibling; // Save next Sibling
		temp->remove();
		minRoot->addSibling(temp);
		temp = nextTemp;
	}

	// Select the left-most sibling of minRoot
	temp = minRoot->leftMostSibling();

	// Remove minRoot and set it to any sibling, if there exists one
	if(temp == minRoot)
	{
		if(minRoot->rightSibling)
			temp = minRoot->rightSibling;
		else
		{
			// Heap is obviously empty
			Node* out = minRoot;
			minRoot->remove();
			minRoot = NULL;
			return out;
		}
	}
	Node* out = minRoot;
	minRoot->remove();
	minRoot = temp;

	// Initialize list of roots	
	for(int i=0; i<100; i++)
		rootListByRank[i] = NULL;
	
	while(temp)
	{
		// Check if key of current vertex is smaller than the key of minRoot
		if(temp->key < minRoot->key)
		{
			minRoot = temp;
		}

		nextTemp = temp->rightSibling;		
		link(temp);
		temp = nextTemp;
	}

	return out;	
}

bool FibonacciHeap::link(Node* root)
{
	// Insert Vertex into root list
	if(rootListByRank[root->rank] == NULL)
	{
		rootListByRank[root->rank] = root;
		return false;
	}
	else
	{
		// Link the two roots
		Node* linkVertex = rootListByRank[root->rank];
		rootListByRank[root->rank] = NULL;
		
		if(root->key < linkVertex->key || root == minRoot)
		{
			linkVertex->remove();
			root->addChild(linkVertex);
			if(rootListByRank[root->rank] != NULL)
				link(root);
			else
				rootListByRank[root->rank] = root;
		}
		else
		{
			root->remove();
			linkVertex->addChild(root);
			if(rootListByRank[linkVertex->rank] != NULL)
				link(linkVertex);
			else
				rootListByRank[linkVertex->rank] = linkVertex;
		}
		return true;
	}
}


void FibonacciHeap::decreaseKey(double delta, Node* vertex)
{
	vertex->key = delta;

	if(vertex->parent != NULL) // The vertex has a parent
	{
		// Remove vertex and add to root list
		vertex->remove();
		minRoot->addSibling(vertex);		
	}
	// Check if key is smaller than the key of minRoot
	if(vertex->key < minRoot->key)
		minRoot = vertex;
}
