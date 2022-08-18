// graph.h <Starter Code>
// Divyam Dhanuka, Fall 2021, UIC
//
// Basic graph class using adjacency matrix representation.  Currently
// limited to a graph with at most 100 vertices.
//
// University of Illinois at Chicago
// CS 251: Fall 2021
// Project #7 - Openstreet Maps
//

#pragma once

#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

using namespace std;

template <typename VertexT, typename WeightT>
class graph {
 private:
  typedef pair<VertexT, WeightT> Edge;
  map<VertexT, set<Edge>> adjList;
  int size;

 public:
  //
  // constructor:
  //
  // Constructs an empty graph where n is the max # of vertices
  // you expect the graph to contain.
  //
  // NOTE: the graph is implemented using an adjacency matrix.
  // If n exceeds the dimensions of this matrix, an exception
  // will be thrown to let you know that this implementation
  // will not suffice.
  //
  graph() { size = 0; }

  //
  // NumVertices
  //
  // Returns the # of vertices currently in the graph.
  //
  int NumVertices() const { return this->adjList.size(); }

  //
  // NumEdges
  //
  // Returns the # of edges currently in the graph.
  //
  int NumEdges() const { return size; }

  //
  // addVertex
  //
  // Adds the vertex v to the graph if there's room, and if so
  // returns true.  If the graph is full, or the vertex already
  // exists in the graph, then false is returned.
  //
  bool addVertex(VertexT v) {
    //
    // is the vertex already in the graph?  If so, we do not
    // insert again otherwise Vertices may fill with duplicates:
    //
    if (this->adjList.count(v) > 0) {
      return false;
    }

    //
    // if we get here, vertex does not exist so insert.
    // Using the operator[] to insert the key with the
    // default empty value into the list
    //
    this->adjList[v];

    return true;
  }

  //
  // addEdge
  //
  // Adds the edge (from, to, weight) to the graph, and returns
  // true.  If the vertices do not exist or for some reason the
  // graph is full, false is returned.
  //
  // NOTE: if the edge already exists, the existing edge weight
  // is overwritten with the new edge weight.
  //
  bool addEdge(VertexT from, VertexT to, WeightT weight) {
    // If either of the vertice do not exist then return false
    if (this->adjList.count(from) == 0 || this->adjList.count(to) == 0)
      return false;

    auto it = adjList.find(from);
    Edge edge(to, weight);

    // Finding if an edge already exists
    for (auto e : it->second) {
      if (e.first == to) {
        // Since edge exists we update the weight for the edge
        Edge edge2(e.first, e.second);
        it->second.erase(edge2);
        size--;
        break;
      }
    }

    // edge does not exist so we insert the new edge into the set for the
    // vertice
    it->second.insert(edge);
    size++;
    return true;
  }

  //
  // getWeight
  //
  // Returns the weight associated with a given edge.  If
  // the edge exists, the weight is returned via the reference
  // parameter and true is returned.  If the edge does not
  // exist, the weight parameter is unchanged and false is
  // returned.
  //
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
    // If either of the vertice do not exist then return false
    if (this->adjList.count(from) == 0 || this->adjList.count(to) == 0)
      return false;

    auto it = adjList.find(from);

    for (auto e : it->second) {
      if (e.first == to) {
        // Since edge exists we update the weight for the edge
        weight = e.second;
        return true;
      }
    }

    // edge does not exist so we return false
    return false;
  }

  //
  // neighbors
  //
  // Returns a set containing the neighbors of v, i.e. all
  // vertices that can be reached from v along one edge.
  // Since a set is returned, the neighbors are returned in
  // sorted order; use foreach to iterate through the set.
  //
  set<VertexT> neighbors(VertexT v) const {
    set<VertexT> S;

    // finding the vertex in the adjacent list
    auto itr = this->adjList.find(v);

    if (itr == adjList.end()) {
      return S;
    }

    for (auto e : itr->second) {
      // looping through all the edges for that vertex and adding them to a set
      S.insert(e.first);
    }

    return S;
  }

  //
  // getVertices
  //
  // Returns a vector containing all the vertices currently in
  // the graph.
  //
  vector<VertexT> getVertices() const {
    vector<VertexT> V;

    for (auto i : this->adjList) {
      // Looping through all the elemnts in the map that only contains the
      // vertices and adding that to a vector.
      V.push_back(i.first);
    }

    return V;
  }

  //
  // dump
  //
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    graph<string,int>  G(26);
  //    ...
  //    G.dump(cout);  // dump to console
  //
  void dump(ostream& output) const {
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;

    output << "**Num vertices: " << this->NumVertices() << endl;
    output << "**Num edges: " << this->NumEdges() << endl;

    output << endl;
    output << "**Vertices:" << endl;
    int counter = 0;
    for (auto j : this->adjList) {
      output << " " << counter << ". " << j.first << endl;
      counter++;
    }

    output << endl;
    output << "**Edges:" << endl;
    int row = 0;
    for (auto i : this->adjList) {
      output << " row " << row << ": ";
      row++;

      for (auto e : i.second) {
        WeightT weight;
        if (!getWeight(i.first, e.first, weight)) {
          output << "F ";
        } else {
          output << "(T," << weight << ") ";
        }
      }
      output << endl;
    }
    output << "**************************************************" << endl;
  }
};
