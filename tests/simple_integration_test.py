#!/usr/bin/env python
"""
Simplified integration tests for panfeed CLI
"""

import os
import sys
import shutil
import subprocess
import filecmp


def run_panfeed_command(args, log_file):
    """Run panfeed command and capture output"""
    cmd = [sys.executable, "panfeed-runner.py"] + args

    # Run command and capture output
    with open(log_file, 'w') as f:
        result = subprocess.run(cmd, cwd=".", stdout=f, stderr=subprocess.STDOUT)

    return result.returncode


def create_absolute_path_files(test_files_dir):
    """Create input files with absolute paths for fileoffiles test"""
    # Check if test_files directory exists
    if not os.path.exists(test_files_dir):
        raise FileNotFoundError(
            f"Test files directory not found: {test_files_dir}\n"
            f"Please extract the test files first by running:\n"
            f"  cd tests && tar -xf test_files.tar.gz"
        )
    
    # Check if required subdirectories exist
    gff_dir = os.path.join(test_files_dir, "gffs")
    fasta_dir = os.path.join(test_files_dir, "fastas")
    
    if not os.path.exists(gff_dir):
        raise FileNotFoundError(
            f"GFF files directory not found: {gff_dir}\n"
            f"Please ensure test_files.tar.gz has been properly extracted"
        )
    
    if not os.path.exists(fasta_dir):
        raise FileNotFoundError(
            f"FASTA files directory not found: {fasta_dir}\n"
            f"Please ensure test_files.tar.gz has been properly extracted"
        )
    
    # List actual files in directories (don't hardcode filenames)
    gff_files = [f for f in os.listdir(gff_dir) if f.endswith('.gff')]
    fasta_files = [f for f in os.listdir(fasta_dir) if f.endswith('.fasta')]
    
    if not gff_files:
        raise FileNotFoundError(
            f"No GFF files found in: {gff_dir}\n"
            f"Please ensure test_files.tar.gz has been properly extracted"
        )
    
    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA files found in: {fasta_dir}\n"
            f"Please ensure test_files.tar.gz has been properly extracted"
        )
    
    # Create absolute path files
    gff_abs_file = os.path.join(test_files_dir, "input_gffs_abs.txt")
    fasta_abs_file = os.path.join(test_files_dir, "input_fastas_abs.txt")
    
    with open(gff_abs_file, 'w') as f:
        for gff_file in gff_files:
            abs_path = os.path.abspath(os.path.join(gff_dir, gff_file))
            f.write(f"{abs_path}\n")
    
    with open(fasta_abs_file, 'w') as f:
        for fasta_file in fasta_files:
            abs_path = os.path.abspath(os.path.join(fasta_dir, fasta_file))
            f.write(f"{abs_path}\n")

def run_integration_tests():
    """Run comprehensive integration tests"""

    # Test configuration
    # Use dynamic path resolution that works from both root and tests directory
    if os.path.exists("test_files"):
        test_files_dir = "test_files"
    else:
        test_files_dir = "tests/test_files"
    test_output_base = os.path.join(test_files_dir, "simple_integration_test_output")

    # Create test output directory
    if os.path.exists(test_output_base):
        shutil.rmtree(test_output_base)
    os.makedirs(test_output_base, exist_ok=True)

    # Create absolute path files for fileoffiles test
    create_absolute_path_files(test_files_dir)

    # Test scenarios
    test_scenarios = [
        {
            'name': 'basic',
            'description': 'Basic run with default parameters',
            'args': [
                '--gff', os.path.join(test_files_dir, 'gffs'),
                '--presence-absence', os.path.join(test_files_dir, 'gene_presence_absence.csv'),
                '--targets', os.path.join(test_files_dir, 'stroi.txt'),
                '--output', os.path.join(test_output_base, 'basic')
            ]
        },
        {
            'name': 'upstream_downstream',
            'description': 'Upstream and downstream regions',
            'args': [
                '--gff', os.path.join(test_files_dir, 'gffs'),
                '--presence-absence', os.path.join(test_files_dir, 'gene_presence_absence.csv'),
                '--targets', os.path.join(test_files_dir, 'stroi.txt'),
                '--upstream', '100',
                '--downstream', '100',
                '--output', os.path.join(test_output_base, 'upstream_downstream')
            ]
        },
        {
            'name': 'noncanonical',
            'description': 'Non-canonical k-mers',
            'args': [
                '--gff', os.path.join(test_files_dir, 'gffs'),
                '--presence-absence', os.path.join(test_files_dir, 'gene_presence_absence.csv'),
                '--targets', os.path.join(test_files_dir, 'stroi.txt'),
                '--non-canonical',
                '--output', os.path.join(test_output_base, 'noncanonical')
            ]
        },
        {
            'name': 'compressed',
            'description': 'Compressed output',
            'args': [
                '--gff', os.path.join(test_files_dir, 'gffs'),
                '--presence-absence', os.path.join(test_files_dir, 'gene_presence_absence.csv'),
                '--targets', os.path.join(test_files_dir, 'stroi.txt'),
                '--compress',
                '--output', os.path.join(test_output_base, 'compressed')
            ]
        },
        {
            'name': 'fileoffiles',
            'description': 'File of files input',
            'args': [
                '--gff', os.path.join(test_files_dir, 'input_gffs_abs.txt'),
                '--fasta', os.path.join(test_files_dir, 'input_fastas_abs.txt'),
                '--presence-absence', os.path.join(test_files_dir, 'gene_presence_absence.csv'),
                '--targets', os.path.join(test_files_dir, 'stroi.txt'),
                '--output', os.path.join(test_output_base, 'fileoffiles')
            ]
        }
    ]

    # Run all test scenarios
    results = []

    for scenario in test_scenarios:
        print(f"\n{'='*60}")
        print(f"Running test: {scenario['name']}")
        print(f"Description: {scenario['description']}")
        print(f"{'='*60}")

        # Create log file path
        log_file = os.path.join(test_output_base, f"{scenario['name']}.log")

        # Clean up output directory if it exists
        output_dir = scenario['args'][-1]
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        # Run the test
        return_code = run_panfeed_command(scenario['args'], log_file)

        # Check if command succeeded
        if return_code != 0:
            print(f"❌ Test {scenario['name']} FAILED with return code {return_code}")
            results.append({
                'name': scenario['name'],
                'status': 'FAILED',
                'return_code': return_code
            })
            continue

        print(f"✅ Test {scenario['name']} completed successfully")
        results.append({
            'name': scenario['name'],
            'status': 'PASSED',
            'return_code': return_code
        })

    # Print summary
    print(f"\n{'='*60}")
    print("INTEGRATION TEST SUMMARY")
    print(f"{'='*60}")

    passed = sum(1 for r in results if r['status'] == 'PASSED')
    failed = sum(1 for r in results if r['status'] == 'FAILED')

    print(f"Total tests: {len(results)}")
    print(f"✅ Passed: {passed}")
    print(f"❌ Failed: {failed}")

    # Print detailed results
    print(f"\nDetailed Results:")
    for result in results:
        status_symbol = '✅' if result['status'] == 'PASSED' else '❌'
        print(f"  {status_symbol} {result['name']}: {result['status']}")

    # Save results to file
    results_file = os.path.join(test_output_base, 'simple_integration_test_results.json')
    import json
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to: {results_file}")

    # Return overall status
    return failed == 0  # True if no failures


if __name__ == "__main__":
    # Change to project root directory
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    success = run_integration_tests()

    if success:
        print("\n🎉 All integration tests completed successfully!")
        sys.exit(0)
    else:
        print("\n💥 Some integration tests failed!")
        sys.exit(1)
