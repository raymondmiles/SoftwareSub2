# Code organisation
├── example_file.fasta
├── fastaTogenbank.py
├── __init__.py
├── OUTPUTS
├── README.md
├── requirements.txt
├── test_script.py
└── workingtree.txt
    
# How to run the script

The module fastaTogenbank produces a genbank file from a fasta file. This script is run through the command:
'''
python3 fastaTogenbank.py -i example_file.fasta -o OUTPUTS/example_file.gb
'''
The prefix -i or --input is a required argument which specifies the fasta file that will be converted. 

The prefix -o or --output allows the user to pass the optional argument that is the path and the name of the genbank file that is being created. This is an optional argument and if it is not specified then genbankOutput will be the default name of the file created within the current directory.


If the optional input -p or --print_stats argument is passed then the following statistics of the fasta file are printed to the command line:
    Mean read length
    Median read length
    Maximum read length
    Minimum read length
    Total number of reads


# Running the unit tests

This command runs the test_script.py script and collects data on code coverage:
'''
coverage run test_script.py
'''
This command generates a summary report in the command line that shows the coverage percentage for each file tested:
'''
coverage report
'''
This command generates a html report with a visual representation of coverage. You can access the report within a directory named htmlcov and view it by opening the index.html file in the web browser:
'''
coverage html
'''

# Help
You can run this command to get more information on running the code
'''
python3 fastaTogenbank.py -h
'''
