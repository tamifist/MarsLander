using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;

class Program
{
    static void Main(string[] args)
    {
        int surfaceN = int.Parse(Console.ReadLine()); // the number of points used to draw the surface of Mars.
        Point[] land = new Point[surfaceN];

        for (int i = 0; i < surfaceN; i++)
        {
            string[] inputs = Console.ReadLine().Split(' ');
            land[i] = new Point(int.Parse(inputs[0]), int.Parse(inputs[1]));
        }
        Surface surface = new Surface(land, surfaceN);

        Shuttle shuttle = new Shuttle(surface);
        int turnNb = 0;
        GeneticAlgorithm geneticAlgorithm = new GeneticAlgorithm(30, 100, 0.9, 0.01, 3, 3, 10000, 100, surface, shuttle);

        // game loop
        while (true)
        {
            string[] inputs = Console.ReadLine().Split(' ');
            int x = int.Parse(inputs[0]);
            int y = int.Parse(inputs[1]);
            int hSpeed = int.Parse(inputs[2]); // the horizontal speed (in m/s), can be negative.
            int vSpeed = int.Parse(inputs[3]); // the vertical speed (in m/s), can be negative.
            int fuel = int.Parse(inputs[4]); // the quantity of remaining fuel in liters.
            int angle = int.Parse(inputs[5]); // the rotation angle in degrees (-90 to 90).
            int power = int.Parse(inputs[6]); // the thrust power (0 to 4).

            // Update the shuttle status;
            if (turnNb > 0)
            {
                shuttle.PerformOneTurn(geneticAlgorithm.ToPlay.TiltAngle, geneticAlgorithm.ToPlay.ThrustPower);
            }
            else
            {
                shuttle.CurrentLocation = new Point(x, y);
                shuttle.HorizontalSpeed = hSpeed; // the horizontal speed (in m/s), can be negative.
                shuttle.VerticalSpeed = vSpeed; // the vertical speed (in m/s), can be negative.
                shuttle.RemainingFuel = fuel; // the quantity of remaining fuel in liters.
                shuttle.CurrentFlightInstruction = new FlightInstruction(angle, power);
                shuttle.FlightStatus = FlightStatus.Flying;
            }

            geneticAlgorithm.Run();

            Console.WriteLine(geneticAlgorithm.ToPlay.TiltAngle + " " + geneticAlgorithm.ToPlay.ThrustPower); // 2 integers: rotate power. rotate is the desired rotation angle (should be 0 for level 1), power is the desired thrust power (0 to 4).
            turnNb++;
        }
    }
}

public class GeneticAlgorithm
{
    private readonly int populationSize;
    private readonly int chromosomeLength;
    private readonly double crossoverRate;
    private readonly double mutationRate;
    private readonly int elitismCount;
    private readonly int tournamentSize;
    private readonly int maxGenerations;
    private readonly int timeout;
    private readonly Surface surface;
    private readonly Shuttle shuttle;

    public GeneticAlgorithm(int populationSize, int chromosomeLength, double crossoverRate,
        double mutationRate, int elitismCount, int tournamentSize, int maxGenerations, int timeout, Surface surface, Shuttle shuttle)
    {
        this.populationSize = populationSize;
        this.chromosomeLength = chromosomeLength;
        this.crossoverRate = crossoverRate;
        this.mutationRate = mutationRate;
        this.elitismCount = elitismCount;
        this.tournamentSize = tournamentSize;
        this.maxGenerations = maxGenerations;
        this.timeout = timeout;
        this.surface = surface;
        this.shuttle = shuttle;
    }

    public FlightInstruction ToPlay { get; set; }

    public void Run()
    {
        Stopwatch stopwatch = new Stopwatch();
        stopwatch.Start();

        Population population = new Population(populationSize, chromosomeLength);

        EvalPopulation(population);

        int currentGeneration = 1;
        while (!IsTerminationConditionMet(currentGeneration, stopwatch))
        {
            //population.Individuals = population.OrderByFitness();

            population = CrossoverPopulation(population);
            population = MutatePopulation(population);

            EvalPopulation(population);

            currentGeneration++;
        }
    }
    
    private Population MutatePopulation(Population population)
    {
        Population newPopulation = new Population(population.Individuals.Length);

        for (int individualIndex = 0; individualIndex < population.Individuals.Length; individualIndex++)
        {
            Individual individual = population.GetFittest(individualIndex);//population.Individuals[individualIndex];

            if (individualIndex >= elitismCount)
            {
                for (int geneIndex = 0; geneIndex < individual.Chromosome.Length; geneIndex++)
                {
                    if (Utils.RandomDouble() < mutationRate)
                    {
                        individual.Chromosome[geneIndex] = Utils.GetRandomFlightInstruction();
                    }
                }
            }
            
            newPopulation.Individuals[individualIndex] = individual;
        }

        return newPopulation;
    }

    //void selectionByTournament()
    //{
    //    int ind1, ind2, winner;
    //    int selectedIndNb = 0;
    //    for (int indIdx = 0; indIdx < SELECTED_NB; indIdx++)
    //    {
    //        ind1 = Utils.RandomInt(0, POPULATION_SIZE);
    //        ind2 = Utils.RandomInt(0, POPULATION_SIZE);
    //        if (this.fitness[ind1] > this.fitness[ind2])
    //        {
    //            winner = ind1;
    //        }
    //        else
    //        {
    //            winner = ind2;
    //        }
    //        Utils.copyContent(this.currentPopulation, winner, this.nextPopulation, selectedIndNb, this.individualEffectiveSize);
    //        selectedIndNb++;
    //    }
    //    // Elitism ! (override 1st prev selected)
    //    int bestInd = 0;
    //    int bestFit = this.fitness[0];
    //    for (int i = 1; i < POPULATION_SIZE; i++)
    //    {
    //        if (this.fitness[i] > bestFit)
    //        {
    //            bestFit = this.fitness[i];
    //            bestInd = i;
    //        }
    //    }

    //    Utils.copyContent(currentPopulation, bestInd, this.nextPopulation, 0, this.individualEffectiveSize);
    //}

    private Population CrossoverPopulation(Population population)
    {
        Population newPopulation = new Population(population.Individuals.Length);
        for (int individualIndex = 0; individualIndex < population.Individuals.Length; individualIndex++)
        {
            Individual parent1 = population.GetFittest(individualIndex);
            if (Utils.RandomDouble() < crossoverRate && individualIndex >= elitismCount)
            {
                Individual offspring = new Individual(parent1.Chromosome.Length);
                Individual parent2 = SelectParent(population);
                int crossoverPoint1 = Utils.RandomInt(0, parent1.Chromosome.Length);
                int crossoverPoint2 = Utils.RandomInt(0, parent1.Chromosome.Length);
                int temp = crossoverPoint1;
                crossoverPoint1 = Math.Min(crossoverPoint1, crossoverPoint2);
                crossoverPoint2 = Math.Max(temp, crossoverPoint2);
                for (int geneIndex = 0; geneIndex < parent1.Chromosome.Length; geneIndex++)
                {
                    if (geneIndex < crossoverPoint1 || geneIndex >= crossoverPoint2)
                    {
                        offspring.Chromosome[geneIndex] = parent1.Chromosome[geneIndex];
                    }
                    else
                    {
                        offspring.Chromosome[geneIndex] = parent2.Chromosome[geneIndex];
                    }
                }

                newPopulation.Individuals[individualIndex] = offspring;
            }
            else
            {
                newPopulation.Individuals[individualIndex] = parent1;
            }
        }


        return newPopulation;
    }
    
    private Individual SelectParent(Population population)
    {
        Population tournament = new Population(population.Shuffle().Take(tournamentSize).ToArray());
        for (int i = 0; i < tournament.Individuals.Length; i++)
        {
            tournament.Individuals[i] = population.Individuals[i];
        }

        return tournament.GetFittest(0);
    }

    private bool IsTerminationConditionMet(int currentGeneration, Stopwatch stopwatch)
    {
        return currentGeneration > maxGenerations || stopwatch.ElapsedMilliseconds > timeout;
    }

    private void EvalPopulation(Population population)
    {
        Parallel.ForEach(population.Individuals, individual =>
        {
            CalcFitness(individual);
        });

        population.PopulationFitness = 0;
        foreach (Individual individual in population.Individuals)
        {
            population.PopulationFitness += individual.Fitness;
        }

        ToPlay = population.GetFittest(0).Chromosome[0];
    }

    private double CalcFitness(Individual individual)
    {
        Shuttle shuttleSim = new Shuttle(shuttle, surface);
        shuttleSim.PerformLanding(individual.Chromosome);

        double fitness;
        if (shuttleSim.FlightStatus == FlightStatus.Landed)
        {
            fitness = 10 * shuttleSim.RemainingFuel;
        }
        else if (shuttleSim.FlightStatus == FlightStatus.Flying)
        {
            fitness = (int)-surface.HorizontalDistanceFromFlatZone(shuttleSim.CurrentLocation);
            if (!surface.CollideWithGround(shuttleSim.CurrentLocation, new Point(surface.flatMidX, surface.flatY + 1)))
            {
                fitness = 50;
            }
        }
        else
        {
            fitness = (int)(Math.Min(-Math.Abs(shuttleSim.VerticalSpeed) + 40, 0)
                        + Math.Min(-Math.Abs(shuttleSim.HorizontalSpeed) + 20, 0)
                        + (-Math.Abs(shuttleSim.CurrentFlightInstruction.TiltAngle))
                        - surface.HorizontalDistanceFromFlatZone(shuttleSim.CurrentLocation));
        }

        individual.Fitness = fitness;
        return individual.Fitness;
    }
}

public class Population
{
    public Population(int populationSize)
    {
        Individuals = new Individual[populationSize];
        PopulationFitness = 0;
    }

    public Population(int populationSize, int chromosomeLength): this(populationSize)
    {
        for (int i = 0; i < populationSize; i++)
        {
            Individuals[i] = new Individual(chromosomeLength);
        }
    }

    public Population(Individual[] individuals)
    {
        Individuals = individuals;
    }

    public Individual[] Individuals { get; set; }

    public double PopulationFitness { get; set; }

    public Individual GetFittest(int offset)
    {
        return OrderByFitness().ElementAt(offset);
    }

    public Individual[] OrderByFitness()
    {
        return Individuals.OrderByDescending(x => x.Fitness).ToArray();
    }

    public Individual[] Shuffle()
    {
        Individual[] individuals = new List<Individual>(Individuals).ToArray();
        for (int i = Individuals.Length - 1; i > 0; i--)
        {
            int index = Utils.RandomInt(0, i + 1);
            Individual temp = individuals[index];
            individuals[index] = individuals[i];
            individuals[i] = temp;
        }

        return individuals;
    }
}

public class Individual
{
    public Individual(FlightInstruction[] chromosome)
    {
        Chromosome = chromosome;
    }

    public FlightInstruction[] Chromosome { get; set; }

    public Individual(int chromosomeLength)
    {
        Chromosome = new FlightInstruction[chromosomeLength];
        for (int i = 0; i < chromosomeLength; i++)
        {
            Chromosome[i] = Utils.GetRandomFlightInstruction();
        }
    }

    public double Fitness { get; set; }

    public override string ToString()
    {
        return string.Format("[ {0} ]", string.Join(" | ", Chromosome.Select(x  => $"{x.TiltAngle} {x.ThrustPower}" )));
    }
}

public class Shuttle
{
    private readonly Surface surface;

    public Shuttle(Point currentLocation, FlightInstruction currentFlightInstruction, 
        double horizontalSpeed, double verticalSpeed, int remainingFuel, FlightStatus flightStatus, Surface surface)
    {
        this.surface = surface;

        CurrentLocation = new Point(currentLocation.X, currentLocation.Y);
        CurrentFlightInstruction = new FlightInstruction(currentFlightInstruction.TiltAngle, currentFlightInstruction.ThrustPower);
        HorizontalSpeed = horizontalSpeed;
        VerticalSpeed = verticalSpeed;
        RemainingFuel = remainingFuel;
        FlightStatus = flightStatus;
    }

    public Shuttle(Shuttle shuttle, Surface surface)
        : this(shuttle.CurrentLocation, shuttle.CurrentFlightInstruction, shuttle.HorizontalSpeed, shuttle.VerticalSpeed, shuttle.RemainingFuel, shuttle.FlightStatus, surface)
    {
    }

    public Shuttle(Surface surface)
    {
        this.surface = surface;
    }

    public Point CurrentLocation { get; set; }

    public FlightInstruction CurrentFlightInstruction { get; set; }

    public double HorizontalSpeed { get; set; }

    public double VerticalSpeed { get; set; }

    public int RemainingFuel { get; set; }

    public FlightStatus FlightStatus { get; set; }

    public void PerformLanding(FlightInstruction[] flightInstructions)
    {
        int instructionIndex = 0;
        while (FlightStatus == FlightStatus.Flying && instructionIndex < flightInstructions.Length)
        {
            PerformOneTurn(flightInstructions[instructionIndex].TiltAngle, flightInstructions[instructionIndex].ThrustPower);
            instructionIndex++;
        }
    }

    public void PerformOneTurn(int rotate, int power)
    {
        if (power > CurrentFlightInstruction.ThrustPower)
        {
            CurrentFlightInstruction.ThrustPower++;
        }
        else if (power < CurrentFlightInstruction.ThrustPower)
        {
            CurrentFlightInstruction.ThrustPower--;
        }
        if (rotate > CurrentFlightInstruction.TiltAngle)
        {
            CurrentFlightInstruction.TiltAngle = Math.Min(rotate, CurrentFlightInstruction.TiltAngle + 15);
        }
        else if (rotate < CurrentFlightInstruction.TiltAngle)
        {
            CurrentFlightInstruction.TiltAngle = Math.Max(rotate, CurrentFlightInstruction.TiltAngle - 15);
        }
        
        // 2 - Compute remaining fuel
        RemainingFuel -= CurrentFlightInstruction.ThrustPower;
        
        // 3 - Compute new position
        Point prevLocation = new Point(CurrentLocation.X, CurrentLocation.Y);
        CurrentLocation.X = CurrentLocation.X + HorizontalSpeed - 0.5 * (Math.Sin(Utils.ToRadians(CurrentFlightInstruction.TiltAngle)) * CurrentFlightInstruction.ThrustPower);

        CurrentLocation.Y = CurrentLocation.Y + VerticalSpeed + 0.5 * (Math.Cos(Utils.ToRadians(CurrentFlightInstruction.TiltAngle)) * CurrentFlightInstruction.ThrustPower - 3.711);
        
        // 4 - Compute new speed
        HorizontalSpeed = HorizontalSpeed - Math.Sin(Utils.ToRadians(CurrentFlightInstruction.TiltAngle)) * CurrentFlightInstruction.ThrustPower;
        VerticalSpeed = VerticalSpeed + Math.Cos(Utils.ToRadians(CurrentFlightInstruction.TiltAngle)) * CurrentFlightInstruction.ThrustPower - 3.7111;

        // 5 - Compute new status
        if (!surface.CollideWithGround(prevLocation, CurrentLocation))
        {
            FlightStatus = FlightStatus.Flying;
        }
        else
        {
            if (Math.Abs(HorizontalSpeed) <= 20 && Math.Abs(VerticalSpeed) <= 40 && CurrentFlightInstruction.TiltAngle == 0 &&
                surface.CollideWithFlatZone(prevLocation, CurrentLocation))
            {
                FlightStatus = FlightStatus.Landed;
            }
            else
            {
                FlightStatus = FlightStatus.Crashed;
            }
        }
    }
}

public class Surface
{

    public int nbPoints;
    public Point[] points;

    public int flatMinX;
    public int flatMaxX;
    public int flatMidX;
    public int flatY;

    public int flatZoneStartingIndex;

    public Surface(Point[] pts, int nbPoints)
    {
        this.points = pts;
        this.nbPoints = nbPoints;
        int prevY = -1;

        for (int i = 0; i < nbPoints; i++)
        {
            if (points[i].Y == prevY)
            {
                flatMinX = (int)points[i - 1].X;
                flatMaxX = (int)points[i].X;
                flatY = (int)points[i].Y;
                flatMidX = (flatMinX + flatMaxX) / 2;
                this.flatZoneStartingIndex = i - 1;
                break;
            }
            prevY = (int)points[i].Y;
        }
    }

    public bool CollideWithGround(Point p1, Point p2)
    {
        if (!IsIn(p2)) return true;
        bool col = false;
        int seg = 0;
        do
        {
            if (Geometry.Intersect(p1, p2, points[seg], points[seg + 1]))
            {
                col = true;
            }
            seg++;
        } while (!col && seg < nbPoints - 1);
        return col;
    }

    public static bool IsIn(Point p)
    {
        return p.X < 7000 && p.X > 0 && p.Y < 3000;
    }

    public double HorizontalDistanceFromFlatZone(Point p)
    {
        double hDist = 0;
        if (p.X < flatMinX)
        {
            hDist = flatMinX - p.X + 50; // TODO cheater ?
        }
        else if (p.X > flatMaxX)
        {
            hDist = p.X - flatMaxX + 50; // TODO cheater ?
        }
        return hDist;
    }

    double DistanceFromFlatZone(Point p)
    {
        double dist, hDist = 0, vDist;
        if (p.X < flatMinX)
        {
            hDist = flatMinX - p.X + 50; // TODO cheater ?
        }
        else if (p.X > flatMaxX)
        {
            hDist = p.X - flatMaxX + 50; // TODO cheater ?
        }
        vDist = p.Y - flatY;
        dist = Math.Sqrt((hDist * hDist) + (vDist * vDist));
        return dist;
    }

    public bool IsAlignWithTheFlatZone(Point p)
    {
        return p.Y < 3000 && p.X > flatMinX + 50 && p.X < flatMaxX - 50;
    }

    public bool CollideWithFlatZone(Point p1, Point p2)
    {
        return Geometry.Intersect(p1, p2, points[flatZoneStartingIndex], points[flatZoneStartingIndex + 1]);
    }
}

public class FlightInstruction
{
    public FlightInstruction(int tiltAngle, int thrustPower)
    {
        TiltAngle = tiltAngle;
        ThrustPower = thrustPower;
    }

    public int TiltAngle { get; set; }
    public int ThrustPower { get; set; }
}

public enum FlightStatus
{
    Flying = 0,
    Crashed,
    Landed,
}

public class Point
{
    public Point(double x, double y)
    {
        X = x;
        Y = y;
    }

    public double X { get; set; }

    public double Y { get; set; }
}

public static class Utils
{
    //private static readonly int[] TiltAngles = new int[] { -90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90 };
    private static readonly int[] TiltAngles = new int[] { -90, 0, 90 };
    //private static readonly int[] ThrustPower = new[] { 0, 1, 2, 3, 4 };
    private static readonly int[] ThrustPower = new[] { 0, 4, 4 };
    private static readonly Random Random = new Random();

    public static double RandomDouble()
    {
        lock (Random)
        {
            return Random.NextDouble();
        }
    }

    public static int RandomInt(int minValue, int maxValue)
    {
        lock (Random)
        {
            return Random.Next(minValue, maxValue);
        }
    }

    public static int GetRandomTiltAngle()
    {
        return TiltAngles[RandomInt(0, TiltAngles.Length)];
    }

    public static int GetRandomThrustPower()
    {
        return ThrustPower[RandomInt(0, ThrustPower.Length)];
    }

    public static FlightInstruction GetRandomFlightInstruction()
    {
        return new FlightInstruction(GetRandomTiltAngle(), GetRandomThrustPower());
    }

    public static double ToRadians(int degree)
    {
        return Math.PI * degree / 180.0;
    }
}

public enum Orientation
{
    Colinear = 0,
    Clockwise,
    Counterclockwise,
}

public static class Geometry
{
    // The main function that returns true if line segment ‘p1q1’
    // and ‘p2q2’ intersect.
    public static bool Intersect(Point p1, Point q1, Point p2, Point q2)
    {
        // Find the four orientations needed for general and
        // special cases
        Orientation o1 = GetOrientation(p1, q1, p2);
        Orientation o2 = GetOrientation(p1, q1, q2);
        Orientation o3 = GetOrientation(p2, q2, p1);
        Orientation o4 = GetOrientation(p2, q2, q1);
        // General case
        if (o1 != o2 && o3 != o4)
            return true;
        // Special Cases
        // p1, q1 and p2 are colinear and p2 lies on segment p1q1
        if (o1 == Orientation.Colinear && OnSegment(p1, p2, q1)) return true;
        // p1, q1 and p2 are colinear and q2 lies on segment p1q1
        if (o2 == Orientation.Colinear && OnSegment(p1, q2, q1)) return true;
        // p2, q2 and p1 are colinear and p1 lies on segment p2q2
        if (o3 == Orientation.Colinear && OnSegment(p2, p1, q2)) return true;
        // p2, q2 and q1 are colinear and q1 lies on segment p2q2
        if (o4 == Orientation.Colinear && OnSegment(p2, q1, q2)) return true;
        return false; // Doesn’t fall in any of the above cases
    }


    // Given three colinear points p, q, r, the function checks if
    // point q lies on line segment ‘pr’
    private static bool OnSegment(Point p, Point q, Point r)
    {
        if (q.X <= Math.Max(p.X, r.X) && q.X >= Math.Min(p.X, r.X) &&
            q.Y <= Math.Max(p.Y, r.Y) && q.Y >= Math.Min(p.Y, r.Y))
        {
            return true;
        }

        return false;
    }

    // To find orientation of ordered triplet (p, q, r).
    // The function returns following values
    // 0 –> p, q and r are colinear
    // 1 –> Clockwise
    // 2 –> Counterclockwise
    private static Orientation GetOrientation(Point p, Point q, Point r)
    {
        // See 10th slides from following link for derivation of the formula
        // http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
        double val = (q.Y - p.Y) * (r.X - q.X) - (q.X - p.X) * (r.Y - q.Y);
        if (val == 0)
        {
            return Orientation.Colinear;
        }

        return val > 0 ? Orientation.Clockwise : Orientation.Counterclockwise;
    }
}
