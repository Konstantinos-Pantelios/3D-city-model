
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <math.h> 
#include <list>
#include <map>
#include "Point.h"
#include "Rows.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <fstream>
#include <nlohmann/json.hpp>
#include <any>
#include <boost/uuid/uuid.hpp>            
#include <boost/uuid/uuid_generators.hpp> 
#include <boost/uuid/uuid_io.hpp>
using json = nlohmann::json;


using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
typedef K::Point_3                                Point_3;
typedef CGAL::Surface_mesh<Point_3>               Surface_mesh;
typedef Polyhedron_3::Vertex_iterator        Vertex_iterator;
typedef Polyhedron_3::Facet_iterator                   Facet_iterator;
typedef Polyhedron_3::Halfedge_around_facet_circulator Halfedge_facet_circulator;


int read_object(std::string  input,   std::vector<Point> &vertices)
{
  std::cout << "Reading file: " << input << std::endl;
  std::ifstream infile(input.c_str(), std::ifstream::in);
  if (!infile)
  {
    std::cerr << "Input file not found.\n";
    return false;
  }
  std::string cursor;
  std::string line = "";
  std::getline(infile, line);
  while (line != "")
  {
  	if(cursor == "x" | cursor == "y" |cursor == "z" |cursor == "" ){
  		std::getline(infile, line);
  	}
    std::istringstream linestream(line);

    linestream >> cursor;
    
    float x,y,z;
    x = std::stof(cursor);
 
		linestream >> cursor;

          y = std::stof(cursor);
linestream >> cursor;
          z = std::stof(cursor);

        
      
      auto p = Point(x,y,z);
      vertices.push_back(p);
     
    std::getline(infile, line);

    
  }
  return 0;
}

float distance(Point a, Point b){
	float c = (pow((pow((a.x - b.x),2.0) + pow((a.y - b.y),2.0)),0.5)); //only horizontal
	return c;
}


void getTreeTops(std::vector<Point> &vertices, std::vector<Point> &treeTops){
	for(int i = 0; i < vertices.size(); i++){
		bool itop = true;
		for(int j = 0; j < vertices.size(); j++){	
			if(i!=j && distance(vertices[i], vertices[j]) < 4){
				if(vertices[i].z < vertices[j].z){
					itop = false;
					break;
				}
			}
		}
		if (itop == true){
			auto p = Point(vertices[i].x, vertices[i].y, vertices[i].z);
			treeTops.emplace_back(p);
		}
	}
}

void assignPointsToTree(std::vector<Point> &vertices, std::vector<Point> treeTops, std::map<int, std::vector<int>> &Trees){

for(int i = 0; i < vertices.size(); i++){
	
	float mindistance = 15.0; 	//starting off at 15 making it the maximum distance form a tree
	int treenumber = 0;

	for(int j = 0; j < treeTops.size(); j++){

		if(distance(vertices[i], treeTops[j]) < mindistance){

			mindistance = distance(vertices[i], treeTops[j]);

			treenumber = j;

		}

	}

	Trees[treenumber].emplace_back(i);
	
}
}

void write_file(std::vector<Point> &vertices, std::vector<Point> treeTops, 
	std::map<int, std::vector<int>> &Trees, std::string output_file, 
	std::map<int, std::vector<Point_3>> &points3) {
        std:: ofstream MyFile(output_file);


        MyFile << "x y z tree " << " \n";
        int t = 0;
        for (int i = 0; i < Trees.size(); i++) {
        	if (Trees[i].size() > 8){
            for (int j = 0; j < Trees[i].size(); j++) {
            	
                MyFile << vertices[Trees[i][j]].x << " "<< vertices[Trees[i][j]].y << " "<< vertices[Trees[i][j]].z << " "<< i << " \n";
            	
            	points3[t].push_back(Point_3 (vertices[Trees[i][j]].x,
            	vertices[Trees[i][j]].y,
            	vertices[Trees[i][j]].z)
            	
            	);

            }
                t++;

            }
            
            }
        MyFile.close();
    }


void buildHull(std::map<int, std::vector<Point_3>> &points3,   std::map<int, Polyhedron_3> &polyTrees){
	for ( int i = 0; i < points3.size(); i++){
		CGAL::convex_hull_3(points3[i].begin(), points3[i].end(), polyTrees[i]);
	}
}
void buildTrunk(std::vector<Point_3> &points,  Polyhedron_3 &poly){
	for ( int i = 0; i < points.size(); i++){
		CGAL::convex_hull_3(points.begin(), points.end(), poly);
	}
}

void writeJSON(std::map<int, std::vector<Point_3>> &points3, std::map<int, Polyhedron_3> &polyTrees, std::string output_file){
	json trees;
	json vertices;
	int vertexcount = 0;
	json final;
	final["type"] = "CityJSON";
	final["version"] = "1.0";

	for (int i=1; i < polyTrees.size(); i++){	
			//for every tree
		json boundaries;
		json trunkboundaries;
		json TopTrees;
		json treetop;
		for (  Facet_iterator k = polyTrees[i].facets_begin(); k != polyTrees[i].facets_end(); ++k) {
        Halfedge_facet_circulator j = k->facet_begin();
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
       	json face1= json::array();
        do {
        	face1.emplace_back(std::distance(polyTrees[i].vertices_begin(), j->vertex()) + vertexcount);
        		

        } while ( ++j != k->facet_begin());
        
        json nestedFace;
        nestedFace.emplace_back(face1);
        treetop.emplace_back(nestedFace);
       
		}
		boundaries.emplace_back(json::array().emplace_back(json::array().emplace_back(treetop)));

			float minx = 999999999999;
			float maxx = -999999999999;
			float miny = 999999999999;
			float maxy = -999999999999;
			float minz = 999999999999;
			float maxz = -999999999999;
		
		for ( Vertex_iterator v = polyTrees[i].vertices_begin(); v != polyTrees[i].vertices_end(); ++v){
		json point1 = json::array();
		point1.emplace_back(v->point().hx());
		point1.emplace_back(v->point().hy());
		point1.emplace_back(v->point().hz());

		vertices.emplace_back(point1);
		vertexcount = vertexcount + 1;
		if (minx > v->point().hx()){
			minx = v->point().hx();
		}
		if (maxx < v->point().hx()){
			maxx = v->point().hx();
		}
		if (miny > v->point().hy()){
			miny = v->point().hy();
		}
		if (maxy < v->point().hy()){
			maxy = v->point().hy();
		}
		if (minz > v->point().hz()){
			minz = v->point().hz();
		}
		if (maxz < v->point().hz()){
			maxz = v->point().hz();
		}
			}
	Point_3 topTrunk((maxx + minx)/2, (maxy + miny)/2, minz - 0.3);
	Point_3 bottomTrunk((maxx + minx)/2, (maxy + miny)/2, minz - 5);
	json treetrunk;
	std::vector<Point_3> pointsTrunk;
	Polyhedron_3 trunk;
	Point_3 topTrunk1(((maxx + minx)/2) -0.3, ((maxy + miny)/2) -0.3 , minz - 0.3);
	Point_3 topTrunk2(((maxx + minx)/2) -0.3, ((maxy + miny)/2) +0.3 , minz - 0.3);
	Point_3 topTrunk3(((maxx + minx)/2) +0.3, ((maxy + miny)/2) -0.3 , minz - 0.3);
	Point_3 topTrunk4(((maxx + minx)/2) +0.3, ((maxy + miny)/2) +0.3 , minz - 0.3);
	Point_3 bottomTrunk1(((maxx + minx)/2) -0.3, ((maxy + miny)/2) -0.3 , 0);
	Point_3 bottomTrunk2(((maxx + minx)/2) +0.3, ((maxy + miny)/2) -0.3 , 0);
	Point_3 bottomTrunk3(((maxx + minx)/2) -0.3, ((maxy + miny)/2) +0.3 , 0);
	Point_3 bottomTrunk4(((maxx + minx)/2) +0.3, ((maxy + miny)/2) +0.3 , 0);


	pointsTrunk.push_back(topTrunk1);
	pointsTrunk.push_back(topTrunk2);
	pointsTrunk.push_back(topTrunk3);
	pointsTrunk.push_back(topTrunk4);
	pointsTrunk.push_back(bottomTrunk1);
	pointsTrunk.push_back( bottomTrunk2);
	pointsTrunk.push_back(bottomTrunk3);
	pointsTrunk.push_back(bottomTrunk4);
		
	buildTrunk(pointsTrunk, trunk);

	

	json TRunks;
	json trunkjson;
		for (  Facet_iterator k = trunk.facets_begin(); k != trunk.facets_end(); ++k) {
        Halfedge_facet_circulator j = k->facet_begin();
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
       	json face1= json::array();
        do {
        	face1.emplace_back(std::distance(trunk.vertices_begin(), j->vertex()) + vertexcount);

        } while ( ++j != k->facet_begin());
        
        json nestedTrunk;
        nestedTrunk.emplace_back(face1);
        trunkjson.emplace_back(nestedTrunk);
        

		}
		
		trunkboundaries.emplace_back(json::array().emplace_back(json::array().emplace_back(trunkjson)));
		
		for ( Vertex_iterator v = trunk.vertices_begin(); v != trunk.vertices_end(); ++v){
		json point1 = json::array();
		point1.emplace_back(v->point().hx());
		point1.emplace_back(v->point().hy());
		point1.emplace_back(v->point().hz());

		vertices.emplace_back(point1);
		vertexcount = vertexcount + 1;
		}

		json btype;
		btype["type"] = "Solid";
		json blod;
		blod["lod"] = 2;
		boost::uuids::uuid uuid = boost::uuids::random_generator()();
		string ustring = boost::uuids::to_string(uuid);
		trees[ustring]["type"] = "SolitaryVegetationObject";
		json bmap {{"type", "Solid"}, {"lod", 2}, {"boundaries", boundaries} };
		json tmap {{"type", "Solid"}, {"lod", 2}, {"boundaries", trunkboundaries} };


		trees[ustring]["geometry"].emplace_back(bmap);
		trees[ustring]["geometry"].emplace_back(tmap);

	
	}
	
	final["vertices"] = vertices;
	final["CityObjects"] = trees;

		std:: ofstream MyFile(output_file);
        MyFile << final.dump(4);
        MyFile.close();

}

int main()
{
	std::vector<Point> vertices;
	std::vector<Point> treeTops;
	std::vector<Point_3> points;
	std::map<int, std::vector<Point_3>> points3;
	std::map<int, Polyhedron_3> polyTrees;
  	Polyhedron_3 poly;
	
	std::map<int, std::vector<int>> Trees;
	const char *file_in = "BR.xyz";
	std:: string  input =  file_in;
  	input = "../" + input;
  	const char *file_out = "Treesfile.txt";
  	const char *file_out2 = "BR.json";

	std:: string  output =  file_out;
  	output = "../" + output;
  	std:: string  output2 =  file_out2;
  	output2 = "../" + output2;
	read_object(input, vertices);
	getTreeTops(vertices, treeTops);

	assignPointsToTree(vertices, treeTops, Trees);


    write_file(vertices, treeTops, Trees, output, points3);
    buildHull(points3, polyTrees);

  	Polyhedron_3 poly1 = polyTrees[8];
  	writeJSON(points3, polyTrees, output2);

}
