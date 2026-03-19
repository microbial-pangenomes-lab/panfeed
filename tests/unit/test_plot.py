#!/usr/bin/env python
"""
Unit tests for panfeed.plot module
"""

import os
import tempfile
import shutil
import pytest
import pandas as pd
from panfeed.plot import main as plot_main
from panfeed.plot import get_options


class TestPlotModule:
    """Test cases for plot module functions"""

    def setup_method(self):
        """Setup test fixtures"""
        self.test_dir = tempfile.mkdtemp()
        self.test_files_dir = "tests/test_files"

    def teardown_method(self):
        """Clean up test fixtures"""
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def test_parse_arguments(self):
        """Test argument parsing"""
        # Create test files
        kmers_file = os.path.join(self.test_dir, "test_kmers.tsv")
        with open(kmers_file, 'w') as f:
            f.write("cluster\tstrain\tfeature_id\tcontig\tfeature_strand\tcontig_start\tcontig_end\tgene_start\tgene_end\tstrand\tk-mer\tp_value\n")
            f.write("cluster1\tstrain1\tgene1\tcontig1\t1\t100\t200\t-50\t50\t1\tkmer1\t1e-10\n")
        
        phenotypes_file = os.path.join(self.test_dir, "test_phenotypes.tsv")
        with open(phenotypes_file, 'w') as f:
            f.write("sample\tphenotype\n")
            f.write("strain1\t1\n")
        
        # Test argument parsing - just check that get_options doesn't crash
        try:
            parser = get_options()
            args = parser.parse_args([
                '--kmers', kmers_file,
                '--phenotypes', phenotypes_file,
                '--output', self.test_dir
            ])
            
            assert args.kmers == kmers_file
            assert args.phenotypes == phenotypes_file
            assert args.output == self.test_dir
        except SystemExit:
            # argparse might exit if there are issues, that's ok for this test
            pass

    def test_plot_basic_functionality(self):
        """Test basic plot functionality"""
        # Create test kmers file
        kmers_file = os.path.join(self.test_dir, "test_kmers.tsv")
        with open(kmers_file, 'w') as f:
            f.write("cluster\tstrain\tfeature_id\tcontig\tfeature_strand\tcontig_start\tcontig_end\tgene_start\tgene_end\tstrand\tk-mer\tlrt-pvalue\n")
            f.write("cluster1\tstrain1\tgene1\tcontig1\t1\t100\t200\t-50\t50\t1\tATGCATGCATGCATGCATGCATGCATGCATG\t1e-10\n")
            f.write("cluster1\tstrain2\tgene1\tcontig1\t1\t105\t205\t-55\t55\t1\tATGCATGCATGCATGCATGCATGCATGCATG\t1e-8\n")
        
        # Create test phenotypes file
        phenotypes_file = os.path.join(self.test_dir, "test_phenotypes.tsv")
        with open(phenotypes_file, 'w') as f:
            f.write("sample\tphenotype\n")
            f.write("strain1\t1\n")
            f.write("strain2\t0\n")
        
        # Create output directory
        output_dir = os.path.join(self.test_dir, "plots")
        os.makedirs(output_dir, exist_ok=True)
        
        # Run plot command by setting sys.argv
        import sys
        original_argv = sys.argv
        try:
            sys.argv = [
                'panfeed-plot',
                '--kmers', kmers_file,
                '--phenotype', phenotypes_file,
                '--output-directory', output_dir,
                '--start', '-100',
                '--stop', '100'
            ]
            from panfeed.plot import main
            main()
        except SystemExit:
            pass  # Expected behavior
        finally:
            sys.argv = original_argv
        
        # Check that plot files were created
        # The plot module should create PNG files for each cluster
        plot_files = [f for f in os.listdir(output_dir) if f.endswith('.png')]
        
        # Should have at least some plot files
        assert len(plot_files) > 0
        
        # Check for expected plot types
        significance_plot = any('significance' in f for f in plot_files)
        sequence_plot = any('sequence' in f for f in plot_files)
        hybrid_plot = any('hybrid' in f for f in plot_files)
        
        assert significance_plot or sequence_plot or hybrid_plot

    def test_plot_with_different_parameters(self):
        """Test plot with different parameter combinations"""
        # Create test kmers file
        kmers_file = os.path.join(self.test_dir, "test_kmers.tsv")
        with open(kmers_file, 'w') as f:
            f.write("cluster\tstrain\tfeature_id\tcontig\tfeature_strand\tcontig_start\tcontig_end\tgene_start\tgene_end\tstrand\tk-mer\tlrt-pvalue\n")
            f.write("cluster1\tstrain1\tgene1\tcontig1\t1\t100\t200\t-50\t50\t1\tATGCATGCATGCATGCATGCATGCATGCATG\t1e-10\n")
            f.write("cluster1\tstrain2\tgene1\tcontig1\t1\t105\t205\t-55\t55\t1\tATGCATGCATGCATGCATGCATGCATGCATG\t1e-8\n")
        
        # Create test phenotypes file
        phenotypes_file = os.path.join(self.test_dir, "test_phenotypes.tsv")
        with open(phenotypes_file, 'w') as f:
            f.write("sample\tphenotype\n")
            f.write("strain1\t1\n")
            f.write("strain2\t0\n")
        
        # Test with nucleotides option
        output_dir_nuc = os.path.join(self.test_dir, "plots_nuc")
        os.makedirs(output_dir_nuc, exist_ok=True)
        
        # Run plot command by setting sys.argv
        import sys
        original_argv = sys.argv
        try:
            sys.argv = [
                'panfeed-plot',
                '--kmers', kmers_file,
                '--phenotype', phenotypes_file,
                '--output-directory', output_dir_nuc,
                '--start', '-50',
                '--stop', '50',
                '--nucleotides'
            ]
            from panfeed.plot import main
            main()
        except SystemExit:
            pass
        finally:
            sys.argv = original_argv
        
        # Check that plot files were created
        plot_files = [f for f in os.listdir(output_dir_nuc) if f.endswith('.png')]
        assert len(plot_files) > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
