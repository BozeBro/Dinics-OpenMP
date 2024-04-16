import java.util.*;

public class GraphGenerator
{
    public static Random random;
    
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

        public void printViz(int src)
        {
            System.out.println(src + " -> " + destination + " [ label = " + capacity + "];");
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

        public void printViz()
        {
            for (Edge e: this.neighbors.values())
            {
                e.printViz(id);
            }
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
        
        void printViz()
        {
            System.out.println("digraph G {");
            for (int i = 0; i < nodes.size(); i++)
            {
                nodes.get(i).printViz();
            }
            System.out.println("}");
        }
    }

    public static void main(String[] args)
    {
        try
        {
            random = new Random(Integer.parseInt(args[0]));

            if (args[1].equals("frontier"))
                frontierEdges(Integer.parseInt(args[2]), Integer.parseInt(args[3]), Integer.parseInt(args[4]), Integer.parseInt(args[5]));
            else
                randomEdges(Integer.parseInt(args[2]), Integer.parseInt(args[3]), Integer.parseInt(args[4]));
        }
        catch (Exception e)
        {
            System.out.println("Give args in format: [<seed> frontier <verts> <edges> <maxOutEdges> <maxCapacity>] or [<seed> random <verts> <edges> <maxCapacity>]");
            e.printStackTrace();
        }
    }

    public static void randomEdges(int verts, int edges, int maxCap)
    {
        Graph g = new Graph(verts);

        for (int i = 0; i < edges; i++)
        {
            while (true)
            {
                int start = (int) (random.nextDouble() * verts);
                int end = (int) (random.nextDouble() * verts);
                int capacity = (int) (random.nextDouble() * (maxCap - 1) + 1);

                if (start == end)
                    continue;

                Edge e = new Edge(end, capacity);
                g.nodes.get(start).addEdge(e);
                break;
            }
        }

        g.print();
    }

    public static void frontierEdges(int verts, int edges, int maxOutEdges, int maxCap)
    {
        Graph g = new Graph(verts);
        LinkedList<Integer> frontier = new LinkedList<>();
        frontier.add(0);

        for (int i = 0; i < edges;)
        {
            int start = frontier.pop();
            int count = (int) (random.nextDouble() * (maxOutEdges + 1));

            if (frontier.isEmpty() && count == 0)
                count++;

            int added = 0;
            while (added < count)
            {
                int end = (int) (random.nextDouble() * verts);
                int capacity = (int) (random.nextDouble() * (maxCap - 1) + 1);

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
