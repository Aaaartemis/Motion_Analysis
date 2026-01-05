#  Motion Analysis from Noisy GPS & Acceleration Data

This project analyses **real-world GPS and accelerometer data** collected from an urban walking trajectory to study **motion patterns, sensor noise, and algorithmic data processing**.  
The work combines **signal processing**, **custom data structures**, and **numerical optimisation** to extract reliable motion information from noisy measurements.

The project was completed as part of a **Data Structures & Algorithms (DSA)** coursework and is implemented in **MATLAB**.

---

## Project Overview

The system:
- Converts GPS data into **local coordinate frames**
- Reconstructs **2D trajectories**
- Computes **distance, velocity, and acceleration**
- Applies **digital filtering and outlier rejection**
- Classifies motion into **behavioural modes** (e.g. strolling, walking, running)
- Uses a **custom linked list data structure** to track motion-mode transitions
- Applies **least-squares optimisation** to predict trajectory segments

The focus is on **robust algorithmic reasoning** under noisy, imperfect real-world data.

---

## Repository Structure

```text
.
├── projectCode.m
│   Main entry point for the project.
│   Loads data, runs preprocessing, trajectory reconstruction,
│   speed classification, and analysis pipeline.
│
├── marshgate.m
│   Processes GPS trajectory data for the Marshgate route.
│   Handles coordinate conversion and trajectory reconstruction.
│
├── bridge.m
│   Processes GPS trajectory data for the Bridge route.
│   Used for comparison between different walking paths.
│
├── accelertionTests.m
│   Performs acceleration analysis, filtering, and noise evaluation
│   on accelerometer data.
│
├── SpeedMode.m
│   Defines motion speed modes (e.g. strolling, walking, running)
│   and contains logic for behavioural classification.
│
├── LinkedListNode.m
│   Custom linked list implementation used to store and process
│   motion-mode transitions over time.
│

```
## How to Run the Project

**Important Data Notice:**  
The **original data files used for the project are not included or have not been re-tested** in this repository.  
The code is provided as-is and is intended to be run with **compatible GPS and accelerometer datasets** of the same structure.  
Users may need to adapt file paths or data-loading sections to work with their own datasets.
---

### Requirements
- **MATLAB** (no additional toolboxes required)
- All project files located in the **same directory**
- Input data files formatted consistently with the expected structure
---

### Run the Main Script

Open MATLAB, navigate to the project directory, and run:

```matlab
projectCode

### Case specific analyses

The files: 
- bridge
- marshagte 

are specific for the data recorded which hasn't been provided
