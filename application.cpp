// application.cpp <Starter Code>
// Divyam Dhanuka, Fall 2021, UIC
//
// University of Illinois at Chicago
// CS 251: Fall 2021
// Project #7 - Openstreet Maps
//
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip> /*setprecision*/
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <stack>
#include <string>
#include <vector>

#include "dist.h"
#include "graph.h"
#include "osm.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

double INF = numeric_limits<double>::max();

class prioritize {
 public:
  bool operator()(const pair<long long, double>& p1,
                  const pair<long long, double>& p2) const {
    return p1.second > p2.second;
  }
};

//
// Implement your creative component application here
// TO DO: add arguments
//
void creative() {}

//
// find
//
// This function is used to find if the given ID is in the given vector
//
bool find(vector<long long> visited, long long ID) {
  for (auto i : visited) {
    if (i == ID) {
      return true;
    }
  }
  return false;
}

//
// searchBuilding
//
// This function is used to find the building given the string
//
bool searchBuilding(string query, vector<BuildingInfo> Buildings,
                    BuildingInfo& building) {
  for (auto i : Buildings) {
    if (query == i.Abbrev) {
      building = i;
      return true;
    }
    if (i.Fullname.find(query) != string::npos) {
      building = i;
      return true;
    }
  }
  return false;
}

//
// nearestBuilding
//
// This function is used to find the nearest building to given coordinatess
//
BuildingInfo nearestBuilding(Coordinates midpoint,
                             vector<BuildingInfo> Buildings,
                             set<long long> unreachable_buildings) {
  double dist_temp;
  stack<pair<BuildingInfo, double>> distances;
  BuildingInfo nearest_building;
  pair<BuildingInfo, double> pair1(nearest_building, INF);
  distances.push(pair1);
  for (auto i : Buildings) {
    // Looping through all the buildings
    if (unreachable_buildings.find(i.Coords.ID) ==
        unreachable_buildings.end()) {
      // if the building is not in the unreachable buildings set then continue
      dist_temp = distBetween2Points(midpoint.Lat, midpoint.Lon, i.Coords.Lat,
                                     i.Coords.Lon);
      if (dist_temp < distances.top().second) {
        // If the distance for this building is less that the lowest one yet
        // then add into distances
        pair<BuildingInfo, double> pair2(i, dist_temp);
        distances.push(pair2);
      }
    }
  }

  // nearest building is the one with the shortest distance which is at teh top
  // of the stack
  nearest_building = distances.top().first;
  return nearest_building;
}

//
// nearestNode
//
// This function is used to find the node nearest to the given Footway
//
Coordinates nearestNode(BuildingInfo b, vector<FootwayInfo> Footways,
                        map<long long, Coordinates> Nodes) {
  double min = INF;
  double dist;
  long long temp_ID;
  Coordinates co, temp_co;
  for (auto i : Footways) {
    // Loop through all the footways
    for (auto j : i.Nodes) {
      // Loop through all the nodes for that footway
      temp_ID = j;
      temp_co = Nodes[temp_ID];
      dist = distBetween2Points(b.Coords.Lat, b.Coords.Lon, temp_co.Lat,
                                temp_co.Lon);
      if (dist < min) {
        // If this distance is less than the lowest one yet then
        // remember the distance and the coordinates of this node.
        min = dist;
        co = temp_co;
      }
    }
  }
  return co;
}

//
// DijkstraShortestPath
//
// This function is used to write the algorithm for Djkstra Shortest Path.
//
void DijkstraShortestPath(long long startV, graph<long long, double> G,
                          map<long long, double>& distances,
                          map<long long, long long>& predecessors) {
  vector<long long> vertices = G.getVertices();
  vector<long long> visited;
  priority_queue<pair<long long, double>, vector<pair<long long, double>>,
                 prioritize>
      unvisitedQueue;
  for (auto i : vertices) {
    // initializing values
    distances[i] = INF;
    predecessors[i] = -1;
    pair<long long, double> pair1(i, INF);
    unvisitedQueue.push(pair1);
  }

  distances[startV] = 0;
  pair<long long, double> pair2(startV, 0);
  unvisitedQueue.push(pair2);
  long long currentV;
  set<long long> neighbors;
  double edgeWeight, altPathDist;
  bool throw_away;

  while (!unvisitedQueue.empty()) {
    currentV = unvisitedQueue.top().first;
    unvisitedQueue.pop();
    // Deque the min from the unvisited Queue

    if (distances[currentV] == INF) {
      break;
    } else if (find(visited, currentV)) {
      continue;
    } else {
      visited.push_back(currentV);  // visit currentV
    }

    neighbors = G.neighbors(currentV);
    // for each vertex that is adjacent to currentV
    for (auto j : neighbors) {
      throw_away = G.getWeight(currentV, j, edgeWeight);
      if (throw_away == false) {  // There shouldnt be such a case
        return;
      }

      // calc the distance for this edge
      altPathDist = distances[currentV] + edgeWeight;

      // if new path weight is less than the one already calculated then
      if (altPathDist < distances[j]) {
        // store the new path weight
        distances[j] = altPathDist;
        // update the predecessor for this vertex to the new vertex
        predecessors[j] = currentV;
        pair<long long, double> pair3(j, altPathDist);
        // push the new path into the unvisited Queue
        unvisitedQueue.push(pair3);
      }
    }
  }
}

//
// getPath
//
// This function is used to get the path to a given vertex given a map of
// predecessors
//
vector<long long> getPath(map<long long, long long> predecessors,
                          long long endVertex) {
  vector<long long> path;
  stack<long long> temp_path;
  long long currV = endVertex;
  while (currV != -1) {
    // loop until there is no vertex predecessor
    temp_path.push(currV);
    // push the current predecessor into the path vector
    currV = predecessors[currV];
  }

  // The path right now is from the destination to the source
  // This loop will create a stack with the path from the source to the
  // destination
  while (!temp_path.empty()) {
    currV = temp_path.top();
    temp_path.pop();
    path.push_back(currV);
  }

  return path;
}

void application(map<long long, Coordinates>& Nodes,
                 vector<FootwayInfo>& Footways, vector<BuildingInfo>& Buildings,
                 graph<long long, double> G);

//
// output
//
// helper function that outputs
//
void output(BuildingInfo nearest_building, Coordinates nearest_co3) {
  cout << "New destination building:" << endl;
  cout << " " << nearest_building.Fullname << endl;
  cout << " (" << nearest_building.Coords.Lat << ", "
       << nearest_building.Coords.Lon << ")" << endl;
  cout << "Nearest destination node:" << endl;
  cout << " " << nearest_co3.ID << endl;
  cout << " (" << nearest_co3.Lat << ", " << nearest_co3.Lon << ")" << endl
       << endl;
}

//
// output2
//
// Helper function that outputs the path taken by person1 and person2 to the
// given destination
//
void output2(vector<long long> path1, vector<long long> path2, double dist1,
             double dist2) {
  cout << "Person 1's distance to dest: " << dist1 << " miles" << endl;
  cout << "Path: " << path1[0];
  for (long unsigned int k = 1; k < path1.size(); ++k) {
    cout << "->" << path1[k];
  }
  cout << endl << endl;
  cout << "Person 2's distance to dest: " << dist2 << " miles" << endl;
  cout << "Path: " << path2[0];
  for (long unsigned int k = 1; k < path2.size(); ++k) {
    cout << "->" << path2[k];
  }
  cout << endl;
}

//
// application2
//
// This is a recursive function used for MS8-MS11
//
void application2(map<long long, Coordinates>& Nodes,
                  vector<FootwayInfo>& Footways,
                  vector<BuildingInfo>& Buildings, graph<long long, double> G,
                  BuildingInfo building1, BuildingInfo building2,
                  set<long long>& unreachable_buildings) {
  Coordinates Co1 = building1.Coords;
  Coordinates Co2 = building2.Coords;
  Coordinates midpoint =
      centerBetween2Points(Co1.Lat, Co1.Lon, Co2.Lat, Co2.Lon);

  BuildingInfo nearest_building =
      nearestBuilding(midpoint, Buildings, unreachable_buildings);

  Coordinates nearest_co1, nearest_co2, nearest_co3;
  nearest_co1 = nearestNode(building1, Footways, Nodes);
  nearest_co2 = nearestNode(building2, Footways, Nodes);
  nearest_co3 = nearestNode(nearest_building, Footways, Nodes);

  output(nearest_building, nearest_co3);

  map<long long, double> distances1, distances2;
  map<long long, long long> predecessors1, predecessors2;

  DijkstraShortestPath(nearest_co1.ID, G, distances1, predecessors1);
  DijkstraShortestPath(nearest_co2.ID, G, distances2, predecessors2);

  vector<long long> path1, path2;
  path1 = getPath(predecessors1, nearest_co3.ID);
  path2 = getPath(predecessors2, nearest_co3.ID);

  // If the destination is unreachable from either sources then find the next
  // closest destination building
  if (distances1[nearest_co3.ID] >= INF || distances2[nearest_co3.ID] >= INF) {
    cout << "At least one person was unable to reach the destination building. "
            "Finding next closest building..."
         << endl
         << endl;
    unreachable_buildings.insert(nearest_building.Coords.ID);
    application2(Nodes, Footways, Buildings, G, building1, building2,
                 unreachable_buildings);
    return;
  }

  output2(path1, path2, distances1[nearest_co3.ID], distances2[nearest_co3.ID]);

  application(Nodes, Footways, Buildings, G);
}

//
// output2
//
// Helper function that is used to output the info
//
void output1(BuildingInfo building1, BuildingInfo building2,
             BuildingInfo nearest_building, Coordinates nearest_co1,
             Coordinates nearest_co2, Coordinates nearest_co3) {
  cout << "Person 1's point:" << endl;
  cout << " " << building1.Fullname << endl;
  cout << " (" << building1.Coords.Lat << ", " << building1.Coords.Lon << ")"
       << endl;
  cout << "Person 2's point:" << endl;
  cout << " " << building2.Fullname << endl;
  cout << " (" << building2.Coords.Lat << ", " << building2.Coords.Lon << ")"
       << endl;
  cout << "Destination Building:" << endl;
  cout << " " << nearest_building.Fullname << endl;
  cout << " (" << nearest_building.Coords.Lat << ", "
       << nearest_building.Coords.Lon << ")" << endl;
  cout << endl;
  cout << "Nearest P1 node:" << endl;
  cout << " " << nearest_co1.ID << endl;
  cout << " (" << nearest_co1.Lat << ", " << nearest_co1.Lon << ")" << endl;
  cout << "Nearest P2 node:" << endl;
  cout << " " << nearest_co2.ID << endl;
  cout << " (" << nearest_co2.Lat << ", " << nearest_co2.Lon << ")" << endl;
  cout << "Nearest destination node:" << endl;
  cout << " " << nearest_co3.ID << endl;
  cout << " (" << nearest_co3.Lat << ", " << nearest_co3.Lon << ")" << endl;
  cout << endl;
}

//
// Implement your standard application here
// TO DO: add a parameter for the graph you make.
//
void application(map<long long, Coordinates>& Nodes,
                 vector<FootwayInfo>& Footways, vector<BuildingInfo>& Buildings,
                 graph<long long, double> G) {
  string person1Building, person2Building;

  cout << endl << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);

  while (person1Building != "#") {
    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);

    BuildingInfo building1, building2;
    bool isBuilding = searchBuilding(person1Building, Buildings, building1);
    if (!isBuilding) {
      cout << "Person 1's building not found" << endl;
      application(Nodes, Footways, Buildings, G);
      return;
    }
    isBuilding = searchBuilding(person2Building, Buildings, building2);
    if (!isBuilding) {
      cout << "Person 2's building not found" << endl;
      application(Nodes, Footways, Buildings, G);
      return;
    }

    Coordinates Co1 = building1.Coords;
    Coordinates Co2 = building2.Coords;
    Coordinates midpoint =
        centerBetween2Points(Co1.Lat, Co1.Lon, Co2.Lat, Co2.Lon);
    set<long long> unreachable_buildings;

    BuildingInfo nearest_building =
        nearestBuilding(midpoint, Buildings, unreachable_buildings);

    Coordinates nearest_co1, nearest_co2, nearest_co3;
    nearest_co1 = nearestNode(building1, Footways, Nodes);
    nearest_co2 = nearestNode(building2, Footways, Nodes);
    nearest_co3 = nearestNode(nearest_building, Footways, Nodes);

    output1(building1, building2, nearest_building, nearest_co1, nearest_co2,
            nearest_co3);

    map<long long, double> distances1, distances2;
    map<long long, long long> predecessors1, predecessors2;

    DijkstraShortestPath(nearest_co1.ID, G, distances1, predecessors1);
    if (distances1[nearest_co2.ID] >= INF) {
      cout << "Sorry, destination unreachable." << endl;
      application(Nodes, Footways, Buildings, G);
      return;
    }
    DijkstraShortestPath(nearest_co2.ID, G, distances2, predecessors2);

    vector<long long> path1, path2;
    path1 = getPath(predecessors1, nearest_co3.ID);
    path2 = getPath(predecessors2, nearest_co3.ID);

    if (distances1[nearest_co3.ID] >= INF ||
        distances2[nearest_co3.ID] >= INF) {
      cout << "At least one person was unable to reach the destination "
              "building. Finding next closest building..."
           << endl << endl;
      unreachable_buildings.insert(nearest_building.Coords.ID);
      application2(Nodes, Footways, Buildings, G, building1, building2,
                   unreachable_buildings);
      return;
    }

    output2(path1, path2, distances1[nearest_co3.ID],
            distances2[nearest_co3.ID]);

    cout << endl << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
  }
}

int main() {
  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates> Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo> Footways;
  // info about each building, in no particular order
  vector<BuildingInfo> Buildings;
  XMLDocument xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;

  graph<long long, double> G;
  for (auto i : Nodes) {
    G.addVertex(i.first);
  }

  for (auto j : Footways) {
    for (long unsigned int k = 0; k < j.Nodes.size() - 1; ++k) {
      Coordinates co1 = Nodes[j.Nodes[k]];
      Coordinates co2 = Nodes[j.Nodes[k + 1]];

      double weight = distBetween2Points(co1.Lat, co1.Lon, co2.Lat, co2.Lon);

      G.addEdge(co1.ID, co2.ID, weight);
      G.addEdge(co2.ID, co1.ID, weight);
    }
  }

  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  //
  // Menu
  //
  string userInput;
  cout << "Enter \"a\" for the standard application or "
       << "\"c\" for the creative component application> ";
  getline(cin, userInput);
  if (userInput == "a") {
    // TO DO: add argument for the graph you make.
    application(Nodes, Footways, Buildings, G);
  } else if (userInput == "c") {
    // TO DO: add arguments
    creative();
  }
  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}
