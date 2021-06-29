//Program to Assemble the Phi-X174 genome using Overlap Graph.

import java.util.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Scanner;
import java.util.StringTokenizer;

public class GenomeAssemblyOverlap {

	// function to get the input reads of the genome.
	public static Vertex[] getReads() {
		FastReader s = new FastReader();
		Vertex[] graph = new Vertex[6];

		for (int i = 0; i < graph.length; i++) {
			graph[i] = new Vertex(i, s.nextLine());
		}
		return graph;
	}

	// class for a Vertex of a Graph.
	static class Vertex {
		int vertexNum; // id of the vertex.
		String read; // read of the vertex.
		Map<Integer, Integer> edges; // Keys are indexes of adjacent vertices, and Values are length of overlap
										// between the two strings (read).
		boolean found; // if found while traversing the graph.

		public Vertex(int vertexNum, String read) {
			this.vertexNum = vertexNum;
			this.read = read;
			this.edges = new HashMap<Integer, Integer>();
			this.found = false;
		}
	}

	// class for reading the input.
	static class FastReader {
		BufferedReader br;
		StringTokenizer st;

		public FastReader() {
			br = new BufferedReader(new InputStreamReader(System.in));
		}

		String next() {
			while (st == null || !st.hasMoreElements()) {
				try {
					st = new StringTokenizer(br.readLine());
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			return st.nextToken();
		}

		int nextInt() {
			return Integer.parseInt(next());
		}

		long nextLong() {
			return Long.parseLong(next());
		}

		double nextDouble() {
			return Double.parseDouble(next());
		}

		String nextLine() {
			String str = "";
			try {
				str = br.readLine();
			} catch (IOException e) {
				e.printStackTrace();
			}
			return str;
		}
	}

	// function to found ssp
	static String str;

	// Utility function to calculate
	// minimum of two numbers
	static int min(int a, int b) {
		return (a < b) ? a : b;
	}

	// Function to calculate maximum
	// overlap in two given strings
	static int findOverlappingPair(String str1, String str2) {

		// max will store maximum
		// overlap i.e maximum
		// length of the matching
		// prefix and suffix
		int max = Integer.MIN_VALUE;
		int len1 = str1.length();
		int len2 = str2.length();

		// check suffix of str1 matches
		// with prefix of str2
		for (int i = 1; i <= min(len1, len2); i++) {

			// compare last i characters
			// in str1 with first i
			// characters in str2
			if (str1.substring(len1 - i).compareTo(str2.substring(0, i)) == 0) {
				if (max < i) {

					// Update max and str
					max = i;
					str = str1 + str2.substring(i);
				}
			}
		}

		// check prefix of str1 matches
		// with suffix of str2
		for (int i = 1; i <= min(len1, len2); i++) {

			// compare first i characters
			// in str1 with last i
			// characters in str2
			if (str1.substring(0, i).compareTo(str2.substring(len2 - i)) == 0) {
				if (max < i) {

					// pdate max and str
					max = i;
					str = str2 + str1.substring(i);
				}
			}
		}

		return max;
	}

	// Function to calculate smallest
	// string that contains
	// each string in the given set as substring.
	static String findShortestSuperstring(ArrayList<String> arr, int len) {

		// run len-1 times to consider every pair
		while (len != 1) {

			// To store maximum overlap
			int max = Integer.MIN_VALUE;

			// To store array index of strings
			// involved in maximum overlap
			int l = 0, r = 0;

			// to store resultant string after
			// maximum overlap
			String resStr = "";

			for (int i = 0; i < len; i++) {
				for (int j = i + 1; j < len; j++) {

					// res will store maximum
					// length of the matching
					// prefix and suffix str is
					// passed by reference and
					// will store the resultant
					// string after maximum
					// overlap of arr[i] and arr[j],
					// if any.
					int res = findOverlappingPair(arr.get(i), arr.get(j));

					// Check for maximum overlap
					if (max < res) {
						max = res;
						resStr = str;
						l = i;
						r = j;
					}
				}
			}

			// Ignore last element in next cycle
			len--;

			// If no overlap,
			// append arr[len] to arr[0]
			if (max == Integer.MIN_VALUE) {
				String temp = arr.get(0);
				temp += arr.get(len);
			}
			else {

				// Copy resultant string
				// to index l
				arr.set(l, resStr);

				// Copy string at last index
				// to index r
				arr.set(r, arr.get(len));
			}
		}
		return arr.get(0);
	}
	// function to make an Overlap graph.
	private static void makeOverlapGraph(Vertex[] graph) {
		for (int i = 0; i < graph.length; i++) {
			for (int j = 0; j < graph.length; j++) {
				if (j == i) {
					continue;
				}
				char[] str1 = graph[i].read.toCharArray();
				char[] str2 = graph[j].read.toCharArray();
				for (int k = 0; k < str2.length; k++) {
					int length = str1.length - k;
					int error = (int) (length * 0.03); // considering an error of only 3% in the overlap between the
														// reads.

					boolean overlap = true;
					int m = 0;
					for (int l = k; l < str1.length; l++) {
						if (error == 0) {
							overlap = false;
							break;
						}
						if (str1[l] != str2[m]) {
							error--;
						}
						m++;
					}

					if (overlap) {
						graph[i].edges.put(j, length);
						break;
					}
				}
			}
		}
	}

	static int last = -1; // store the index of the last explored vertex.

	// function to find hamiltonian path in the graph using a greedy approach. i.e
	// choosing the next vertex with max overlap length.
	private static String findHamiltonianPath(Vertex[] graph) {
		String genome = "";
		int random = 0;
		int first = random;
		ArrayList<String> path = new ArrayList<String>();

		genome = genome + graph[random].read;
		genome = explore(graph, random, genome);
		path.add(genome);

		for (int i = 0; i < graph.length; i++) {
			if (!graph[i].found) {
				path.add(graph[i].read);
				genome = genome + "->" + graph[i].read;
				last = i;
			}
		}

		if (graph[last].edges.containsKey(first)) {
			int length = graph[last].edges.get(first);
			return genome.substring(length);
		}
		System.out.println(findShortestSuperstring(path, path.size()));
		return genome;
	}

	// function to explore the given vertex index.
	private static String explore(Vertex[] graph, int random, String genome) {
		graph[random].found = true;
		int max = -1;
		int index = -1;
		for (int key : graph[random].edges.keySet()) {
			int temp = graph[random].edges.get(key);
			if (graph[key].found == false && max < temp) {
				max = temp;
				index = key;
			}
		}

		if (index == -1) {
			return genome;
		}

		genome = genome + graph[index].read.substring(max);
		last = index;
		return explore(graph, index, genome);
	}

	// main function to run the program.
	public static void main(String[] args) {
		long startime = System.nanoTime();
		Vertex[] graph = getReads(); // get the reads.
		makeOverlapGraph(graph); // make overlap graph.
		String genome = findHamiltonianPath(graph); // find hamiltonian path.
		System.out.println(genome); // print the genome.
		System.out.println("Running time = " + (System.nanoTime() - startime) + "ns");
	}
}