# OH deceleration and trapping simulation

Simulation of deceleration and trapping of OH molecules in the HyT (hybrid trap) experiment in Willitsch group, University of Basel.

## Quick start

### Run locally
1. Define parameters in `src/Parameters.m`
2. Run `main.m` in Matlab

### Submit jobs to a remote server
1. Make a copy this repository on the server
2. Modify `src/Parameters.m` and `main.m` if needed
3. Create a shell script similar to 'scripts/submit_job.sh'
4. Run `sbatch scripts/submit_job.sh` at the root directory
5. Check `logs/` for the job status

## Project Structure

The directory structure looks like this:

```
├── data <- Data required
│ ├── acc <- Deceleration and trapping fields
│ └── sequences <- Decelerator switching sequences
│
├── logs <- Logs generated
│
├── scripts <- Shell scripts
│
├── src <- Source code
│ ├── Deceleration.m <- Deceleration class
│ ├── Beam.m <- Initial molecular beam class
│ ├── Parameters.m <- Define all relevant parameters
│ └── Trapping.py <- Trapping class
│
├── main.m <- Main entry of the program
├── README.md
├── .github <- Github Actions workflows
└── .gitignore <- List of files ignored by git

```

## Features
- Object-oriented programming (OOP) in Matlab
- Support for "parfor" of Matlab to run parallelly
- Support for optional name-value pair arguments in `src/Parameters`, which can be useful for, e.g., optimization tasks. Usage:
`params = Parameters('Trap_coil_current', 88, ...)`