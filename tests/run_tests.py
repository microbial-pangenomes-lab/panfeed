#!/usr/bin/env python
"""
Test runner for panfeed unit tests
"""

import sys
import os
import subprocess
import pytest


def run_unit_tests():
    """Run unit tests"""
    print("Running unit tests...")
    
    # Change to tests directory
    os.chdir("tests")
    
    # Run pytest
    result = pytest.main([
        "unit/",
        "-v",
        "--tb=short",
        "--color=yes"
    ])
    
    return result


def run_integration_tests():
    """Run integration tests"""
    print("Running integration tests...")
    
    # Change to tests directory
    os.chdir("tests")
    
    # Run the existing shell script integration tests
    try:
        result = subprocess.run(["./unit_test.sh"], shell=True, check=True)
        return 0
    except subprocess.CalledProcessError as e:
        print(f"Integration tests failed with exit code {e.returncode}")
        return e.returncode


def main():
    """Main test runner"""
    print("=" * 60)
    print("Panfeed Test Suite")
    print("=" * 60)
    
    # Run unit tests
    unit_result = run_unit_tests()
    
    if unit_result != 0:
        print(f"\nUnit tests failed with exit code {unit_result}")
        sys.exit(unit_result)
    
    # Run integration tests
    integration_result = run_integration_tests()
    
    if integration_result != 0:
        print(f"\nIntegration tests failed with exit code {integration_result}")
        sys.exit(integration_result)
    
    print("\n" + "=" * 60)
    print("All tests passed successfully!")
    print("=" * 60)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
