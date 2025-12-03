# Start Here: Workflow for Developing ML Parameterizations in CESM-MLe

This guide outlines the recommended workflow for creating and integrating a machine-learning–based parameterization within the Community Earth System Model (CESM-MLe). The steps below summarize the end-to-end process leveraging the tools and examples here.

## 1. Identify the Target Bias or Process

Define the model bias or physical process requiring improved representation. Determine how this process is currently formulated in CESM (see the relevant [ESCOMP](https://github.com/ESCOMP) repositories).

## 2. Develop the ML Parameterization

Design and train the ML parameterization, preferably using PyTorch. Ensure that model inputs and outputs map directly to variables exchanged through the corresponding CESM subroutines. Export the trained model to TorchScript to enable inference through FTorch within the CESM runtime environment.

```# Export to TorchScript format (required for FTorch)
model_scripted = torch.jit.script(model)
model_scripted.save("example_model.pt")
```

## 3. Perform Functional Unit Testing

Validate the TorchScript parameterization using a functional unit test (see the CESM-MLe functional testing documentation). The goal is to confirm correct model loading, compatibility, and successful execution through the FTorch interface.

## 4. Integrate Into CESM

Once functional testing is complete, incorporate the parameterization into the appropriate CESM component. An example can be found here. 

## 5. Verify and Iterate

Run controlled test simulations to evaluate scientific performance, numerical stability, and computational cost. Refine the ML parameterization or retrain as needed based on diagnostic outcomes. Consider recalibrating since other model processes may have been compensating for biases in the process you've replaced. 

## 6. Reporting Issues

If challenges arise during any stage of development—ranging from model preparation to functional testing or CESM integration—please file an issue on this GitHub repository here. When doing so, describe the problem, provide minimal reproducible examples if possible, and include relevant logs or error messages to facilitate troubleshooting.
