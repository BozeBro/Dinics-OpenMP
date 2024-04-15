import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;

public class GraphGenerator
{
    static class Edge
    {
        int destination;
        int capacity;

        public Edge(int destination, int capacity)
        {
            this.destination = destination;
            this.capacity = capacity;
        }

        public void print()
        {
            System.out.print(destination + " " + capacity + " ");
        }
    }

    static class Node
    {
        int id;
        HashMap<Integer, Edge> neighbors = new HashMap<>();

        public Node(int id)
        {
            this.id = id;
        }

        void addEdge(Edge e)
        {
            this.neighbors.put(e.destination, e);
        }

        public void print()
        {
            System.out.print(this.neighbors.size() + " ");
            for (Edge e: this.neighbors.values())
            {
                e.print();
            }
            System.out.println();
        }
    }

    static class Graph
    {
        ArrayList<Node> nodes = new ArrayList<>();

        public Graph(int count)
        {
            for (int i = 0; i < count; i++)
            {
                nodes.add(new Node(i));
            }
        }

        void print()
        {
            System.out.println(nodes.size());
            for (int i = 0; i < nodes.size(); i++)
            {
                nodes.get(i).print();
            }
        }
    }

    public static void main(String[] args)
    {
        frontierEdges(10000, 1000000,  2, 100);
    }

    public static void randomEdges(int verts, int edges, int maxCap)
    {
        Graph g = new Graph(verts);

        for (int i = 0; i < edges; i++)
        {
            while (true)
            {
                int start = (int) (Math.random() * verts);
                int end = (int) (Math.random() * verts);
                int capacity = (int) (Math.random() * (maxCap - 1) + 1);

                if (start == end)
                    continue;

                Edge e = new Edge(end, capacity);
                g.nodes.get(start).addEdge(e);
                break;
            }
        }
    }

    public static void frontierEdges(int verts, int edges, int maxOutEdges, int maxCap)
    {
        Graph g = new Graph(verts);
        LinkedList<Integer> frontier = new LinkedList<>();
        frontier.add(0);

        for (int i = 0; i < edges;)
        {
            int start = frontier.pop();
            int count = (int) (Math.random() * (maxOutEdges + 1));

            if (frontier.isEmpty() && count == 0)
                count++;

            int added = 0;
            while (added < count)
            {
                int end = (int) (Math.random() * verts);
                int capacity = (int) (Math.random() * (maxCap - 1) + 1);

                if (start == end)
                    continue;

                frontier.push(end);
                Edge e = new Edge(end, capacity);
                g.nodes.get(start).addEdge(e);
                added++;
                i++;
            }
        }

        g.print();
    }
}
