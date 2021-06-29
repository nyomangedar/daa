//Program to assemble Phi-X174 Genome from its K-Mer Composition.

import java.util.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Scanner;
import java.util.StringTokenizer;

public class GenomeAssemblyDeBruijn {
	// class for reading the k-mers of the genome.
	static class FastReader {
		BufferedReader br;
		StringTokenizer st;

		public FastReader() {
			br = new BufferedReader(new InputStreamReader(System.in));
		}

		String nextLine() {
			String temp = "";
			try {
				temp = br.readLine();
			} catch (IOException e) {
				e.printStackTrace();
			}
			return temp;
		}
	}

	// class for a vertex of the graph.
	static class Vertex {
		int vertexNum; // if of the vertex.
		String str; // string of the vertex.
		List<Integer> edgeList; // list of adjacent vertices.

		public Vertex(int vertexNum, String str, ArrayList<Integer> list) {
			this.vertexNum = vertexNum;
			this.str = str;
			this.edgeList = list;
		}
	}

	// class for a edge of the graph.
	static class Edge {
		int from; // start vertex of the edge.
		int to; // end vertex of the edge.
		boolean used; // check if edge is used while finding eulerian circuit in the graph.

		public Edge(int from, int to) {
			this.from = from;
			this.to = to;
			this.used = false;
		}
	}

	// list of all the edges in the graph.
	static List<Edge> edges;

	// function to read the k-mers as input.
	public static Vertex[] reader() {
		FastReader fr = new FastReader();

		edges = new ArrayList<Edge>();

		HashMap<String, Integer> idmap = new HashMap<String, Integer>(); // map for only unique strings.
		HashMap<String, ArrayList<Integer>> edgesmap = new HashMap<String, ArrayList<Integer>>(); // map for storing the
																									// ids of other
																									// strings which
																									// have an overlap
																									// with the string.

		int id = 0; // to give unique id to the strings.
		String[] testcase = { "tca", "caa", "aaa", "aac", "acg", "cgt", "gtc", "tcg", "cgc", "gcc", "ccg", "cgt", "gtt", "ttt", "ttg", "tgg", "ggc", "gct", "ctg", "tgc", "gcc", "ccc", "ccc", "cca", "cat", "atc", "tct", "ctg", "tgg", "ggc", "gct", "ctt", "ttc", "tcc", "ccc", "ccg", "cga", "gag", "agc", "gca", "cat", "atg",
				"tgg", "ggg", "ggc", "gcc", "ccc", "ccg", "cgc", "gcc", "ccg", "cgt", "gtg", "tgg", "ggg", "ggc", "gcc", "ccc", "cca", "cac", "act", "ctt", "ttc", "tct", "cta", "taa", "aac", "act", "ctc", "tct", "ctg", "tgc", "gcc", "ccg", "cgc", "gcg", "cgg", "ggc", "gct", "ctg", "tgc", "gcg", "cgg", "ggc", "gct", "cta", "tac", "act", "cta", "taa", "aag", "agt", "gtt", "ttg", "tga", "gag", "agt", "gtt", "tta", "taa", };

		// read the k-mer;
		for (int i = 0; i < testcase.length; i++) {
			String str = testcase[i]; // get the kmer.
			String a = str.substring(0, str.length() - 1); // split into two parts. (0-k-1) and (1-k) strings.
			String b = str.substring(1); // 1-k string (a and b are going to be the values of 'str' in the 'vertex'
											// class.)

			if (!idmap.containsKey(a)) {
				idmap.put(a, id);
				id++;
				edgesmap.put(a, new ArrayList<Integer>());
			}

			if (!idmap.containsKey(b)) {
				idmap.put(b, id);
				id++;
				edgesmap.put(b, new ArrayList<Integer>());
			}

			boolean overlap = isOverlap(a, b); // check for overlap. if yes then add id of string b into edgelist of a.
			if (overlap) {
				Edge edge = new Edge(idmap.get(a), idmap.get(b));
				edgesmap.get(a).add(edges.size());
				edges.add(edge);
			}
		}

		Vertex[] graph = maptoArray(idmap, edgesmap); // from the map create the graph. Each vertex's string is either
														// (0-k-1) string or (1-k) string i.e either a or b.
		return graph;
	}

	// function to create the graph from map.
	private static Vertex[] maptoArray(HashMap<String, Integer> idmap, HashMap<String, ArrayList<Integer>> edgesmap) {
		Vertex[] graph = new Vertex[idmap.size()];

		for (String temp : edgesmap.keySet()) {
			graph[idmap.get(temp)] = new Vertex(idmap.get(temp), temp, edgesmap.get(temp));
		}

		return graph;
	}

	// function to check the overlap between the strings.
	private static boolean isOverlap(String a, String b) {
		int j = 0;
		for (int i = 1; i < a.length(); i++) {
			if (a.charAt(i) != b.charAt(j)) {
				return false;
			}
			j++;
		}
		return true;
	}

	// function to find an eulerian cycle in the graph.
	public static String findCycle(Vertex[] graph) {
		List<Integer> cycle = new ArrayList<Integer>();
		explore(graph, 0, cycle);

		// construct the genome from the eulerian cycle.
		String genome = graph[cycle.get(cycle.size() - 1)].str;
		for (int i = cycle.size() - 2; i > -1; i--) {
			String temp = graph[cycle.get(i)].str;
			genome = genome + temp.charAt(temp.length() - 1);

		}

		// printed from 9 because k-mer size is 10 and since it is circular genome so
		// last 9 chars and first 9 chars are same.
		return genome.substring(0);
	}

	// function finds the eulerian cycle.
	private static void explore(Vertex[] graph, int vertex, List<Integer> cycle) {
		List<Integer> edgeList = graph[vertex].edgeList;
		for (int i = 0; i < edgeList.size(); i++) {
			Edge edge = edges.get(edgeList.get(i));
			if (!edge.used) {
				edge.used = true;
				explore(graph, edge.to, cycle);
			}
		}
		cycle.add(vertex);
	}

	// function run. Used it so as to avoid stack overflow problem by using threads.
	public String run() throws IOException {
		Vertex[] graph = reader();
		return findCycle(graph);
	}

	// main function to run the program.
	public static void main(String[] args) throws IOException {
		String genome = "";
			long starttime = System.nanoTime();
			genome = new GenomeAssemblyDeBruijn().run();
			long elapsedtime = System.nanoTime() - starttime;
			System.out.println("Running time = " + elapsedtime + "ns");
		System.out.println("Genome result: " + genome);
	}
}