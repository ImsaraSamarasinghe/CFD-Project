import argparse
from cfd_simulation.model import Model

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Run CFD Simulation.')
    parser.add_argument('--lx', type=float, default=1.0, help='Length in the x-direction')
    parser.add_argument('--ly', type=float, default=1.0, help='Length in the y-direction')
    parser.add_argument('--nx', type=int, default=100, help='Number of grid points in the x-direction')
    parser.add_argument('--ny', type=int, default=100, help='Number of grid points in the y-direction')
    parser.add_argument('--nt', type=int, default=500, help='Number of time steps')
    
    args = parser.parse_args()

    run = Model(args.lx, args.ly, args.nx, args.ny, args.nt)
    run.CreateGrid()
    run.cavity_flow()
    run.plot()

if __name__ == "__main__":
    main()
