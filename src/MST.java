import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class MST {
    public static double totalEdgeWeight(Graph g) {
        double totalWeight = 0;

        for (int i = 0; i < g.numVertices(); i++) {
            // instantiate a new neighbours array with the outneighbours indexes in graph g
            int[] neighbours = g.outNeighbours(i);
            for (int neighbour : neighbours) {
                // ensure each edge is counted only once
                if (i < neighbour) {
                    // calculate the total weight for each neighbour if i is less than neighbour
                    totalWeight = totalWeight + g.weight(i, neighbour);
                }
            }
        }

        return totalWeight;
    }



    public static Graph getRandomGraph(int n, double p) {
        // instantiate a new matrix graph that is not directed
        MatrixGraph graph = new MatrixGraph(n, false);

        for (int i = 0; i < n; ++i) {
            for (int h = i + 1; h < n; ++h) {
                // if random is less than probability p
                if (Math.random() <= p) {
                    double randomWeight = Math.random(); // generate a random weight in the range 0-1
                    // use the graph edge function from MatrixGraph to add index of i and h + random weight to graph
                    graph.addEdge(i, h, randomWeight);
                    // do not add the reverse edge as the graph is undirected
                }
            }
        }

        return graph;
    }




    static class ComponentTracker {
        private int[] representatives;
        private double[] vertexValues;

        public ComponentTracker(int numVertices) {
            representatives = new int[numVertices];
            vertexValues = new double[numVertices];
            // initialise the components with generalised values
            initialiseComponents();
        }

        private void initialiseComponents() {
            for (int i = 0; i < representatives.length; i++) {
                representatives[i] = i;  // each representative is initialised as its own representative
                vertexValues[i] = 0.0; // vertex values are then initialised as 0.0
            }
        }

        int representativeOf(int x) {
            //returns representative of x
            return representatives[x];
        }

        boolean inSameComponent(int x, int y) {
            // returns true if the vertices x and y are in the same component based on representative
            return representativeOf(x) == representativeOf(y);
        }

        void mergeComponents(int x, int y) {
            int componentX = representatives[x];
            int componentY = representatives[y];

            // choose a representative of component x and component y and make a new representative
            int newRepresentative = chooseRepresentative(componentX, componentY);
            // using the new representative update the component with value of componentX and componentY values with new representative from either x or y
            updateComponent(componentX, newRepresentative);
            updateComponent(componentY, newRepresentative);
            // update the vertex value for the new representative based on the getVertexValue function
            vertexValues[newRepresentative] = getVertexValue(newRepresentative);
        }


        int chooseRepresentative(int componentX, int componentY) {
            //if component x is less than component y then  return component x, otherwise return component y
            if (componentX < componentY) {
                return componentX;
            }


            return componentY;
        }

        void updateComponent(int component, int newRepresentative) {
            // for the length of representatives, and if the representative is in the same component then set the value of  that index to the new representative
            // and make vertex value at that index set the getVertexValue index
            for (int i = 0; i < representatives.length; i++) {
                if (representatives[i] == component) {
                    representatives[i] = newRepresentative;
                    vertexValues[i] = getVertexValue(i);
                }
            }

        }

        double getVertexValue(int vertex) {
            // get the vertex value from the vertex array
            return vertexValues[vertex];

        }


        public static Graph minimumSpanningTree(Graph g) {
            int numberVertices = g.numVertices();

            // create an empty graph t with the same number of vertices as g
            Graph T = new MatrixGraph(numberVertices, g.isDirected());
            ComponentTracker componentTracker = new ComponentTracker(numberVertices);
            // to ensure that the mst graph keeps repeating until no more new edges are added
            while (addedgestoMST(g, T, componentTracker)) {

            }
            return T;
        }


        private static boolean addedgestoMST(Graph g, Graph T, ComponentTracker componentTracker) {
            // track whether an edge was added to the mst graph
            boolean edgeadded = false;

            // using a map to store edge information key: "fromVertex-toVertex", value: weight
            Map<String, Double> edges = new HashMap<>();

            // iterate over all vertices in the graph
            for (int fromVertex = 0; fromVertex < g.numVertices(); fromVertex++) {
                int[] neighbours = g.outNeighbours(fromVertex);
                for (int toVertex : neighbours) {
                    // check if the vertices are not in the same connected component
                    if (!componentTracker.inSameComponent(fromVertex, toVertex)) {
                        // if not in the same component then add the edge to the map with its weight
                        double weight = g.weight(fromVertex, toVertex);
                        edges.put(fromVertex + "-" + toVertex, weight);
                    }
                }
            }

            // sort the edges based on weight in ascending order
            ArrayList<Map.Entry<String, Double>> sortedEdges = new ArrayList<>(edges.entrySet());
            sortedEdges.sort(Map.Entry.comparingByValue());

            // iteration over the sorted edges and add them to the mst graph
            for (Map.Entry<String, Double> entry : sortedEdges) {
                // split the key to get the fromVertex and toVertex value
                String[] vertices = entry.getKey().split("-");
                int fromVertex = Integer.parseInt(vertices[0]);
                int toVertex = Integer.parseInt(vertices[1]);

                // check if the vertices are still not in the same connected component
                if (!componentTracker.inSameComponent(fromVertex, toVertex)) {
                    // add the edge to the mst then merge the connected components and set the flag to true
                    T.addEdge(fromVertex, toVertex, entry.getValue());
                    // merge the connected components of its end points
                    componentTracker.mergeComponents(fromVertex, toVertex);
                    edgeadded = true;
                }
            }

            return edgeadded;
        }

        //  use the values from GraphOfEssex and getGraph function
        public static Graph createEssexGraph() {
            return GraphOfEssex.getGraph();
        }

        public static void main(String[] args) {
            // test case 5.1
            Graph essexGraph = createEssexGraph();
            double essexTotalWeight = totalEdgeWeight(essexGraph);

            Graph essexMST = ComponentTracker.minimumSpanningTree(essexGraph);
            double essexMSTWeight = totalEdgeWeight(essexMST);


            // test case 5.2
            Graph randomGraph = getRandomGraph(100, 0.4);
            double randomGraphTotalWeight = totalEdgeWeight(randomGraph);


            // test case 5.3
            int numGraphs = 100;
            int numVertices = 100;
            double edgeProbability = 0.4;
            double totalMSTWeight = 0;
            double averageMSTWeight = 0;

            // while i is less than numGraphs, add numVertices and edgeProbability to the random graph then put that graph
            // into mst to compute the minimum spanning tree
            // then calculate the total edge weight
            for (int i = 0; i < numGraphs; i++) {
                Graph randomGraphInstance = getRandomGraph(numVertices, edgeProbability);
                Graph randomGraphMST = minimumSpanningTree(randomGraphInstance);

                totalMSTWeight = totalMSTWeight + totalEdgeWeight(randomGraphMST);
            }

            averageMSTWeight = totalMSTWeight / numGraphs;


            System.out.println("total edge weight of essex graph :   " + essexTotalWeight);

            System.out.println("total weight of minimum spanning tree for essex graph : " + essexMSTWeight);

            System.out.println("total edge weight of random graph:  " + randomGraphTotalWeight);

            System.out.println("average total weight of minimum spanning trees for random graphs: " + averageMSTWeight);
        }
    }
}