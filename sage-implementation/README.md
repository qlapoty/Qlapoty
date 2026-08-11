# Qlapoty Sage Code

Sage code accompanying the paper "Qlapoty: Improved analysis and efficiency for quaternionic ideal to isogeny transformation".

All code in this repository is written in `Python` and [`SageMath`](https://www.sagemath.org/), and is made publicly available under the MIT License.
Dependencies
- Python 3.11.1
- SageMath version 10.7


## Code Structure

- `qlapoty.py`: Implements the main qlapoty algorithm, including a comparative analysis with the original Qlapoti algorithm, taken from the original [Qlapoti repo](https://github.com/KULeuven-COSIC/Qlapoti).
- `helpers`: Contain helper functions used for the Qlapoti algorithm. They are taken from the original [Qlapoti repo](https://github.com/KULeuven-COSIC/Qlapoti).  
- `tests_heuristics`: Folder containing jupyter notebooks used to test and experimently verify the claims and assuption of the paper. 

## License

This project is licensed under the terms of the MIT License.