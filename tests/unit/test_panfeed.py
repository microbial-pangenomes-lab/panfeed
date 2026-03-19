#!/usr/bin/env python
"""
Unit tests for panfeed.panfeed module
"""

import os
import tempfile
import shutil
import pytest
import numpy as np
import pandas as pd
from io import StringIO
from panfeed.panfeed import (
    init_presabs_vector,
    cluster_cutter,
    pattern_hasher,
    write_headers
)
from panfeed.input import parse_gff, create_faidx, iter_gene_clusters
from panfeed.classes import Seqinfo


class TestPanfeedModule:
    """Test cases for panfeed module functions"""

    def setup_method(self):
        """Setup test fixtures"""
        self.test_dir = tempfile.mkdtemp()
        self.test_files_dir = "tests/test_files"

    def teardown_method(self):
        """Clean up test fixtures"""
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def test_init_presabs_vector(self):
        """Test init_presabs_vector function"""
        n_strains = 5
        clusterpresab = np.array([1, 1, 0, 1, 0])
        
        # Test without missing_nan - should create zeros, not copy clusterpresab
        result = init_presabs_vector(n_strains, clusterpresab, missing_nan=False)
        expected = np.zeros(n_strains, dtype=np.float64)
        np.testing.assert_array_equal(result, expected)
        
        # Test with missing_nan - should create zeros but set missing (0) to NaN
        result_nan = init_presabs_vector(n_strains, clusterpresab, missing_nan=True)
        expected_nan = np.zeros(n_strains, dtype=np.float64)
        expected_nan[clusterpresab == 0] = np.nan
        np.testing.assert_array_equal(result_nan, expected_nan)

    def test_cluster_cutter_basic(self):
        """Test cluster_cutter with basic parameters"""
        # Prepare test data
        presence_absence = os.path.join(self.test_files_dir, "gene_presence_absence.csv")
        genepres = pd.read_csv(presence_absence, sep=",", index_col=0, low_memory=False)
        
        # Get first cluster
        first_cluster = genepres.iloc[0]
        idx = first_cluster.name
        clusterpresab = first_cluster.values
        
        # Prepare genome data for this cluster
        gff_dir = os.path.join(self.test_files_dir, "gffs")
        fasta_dir = os.path.join(self.test_files_dir, "fastas")
        genome_data = {}
        for gff_file in os.listdir(gff_dir):
            if gff_file.endswith(".gff"):
                strain = gff_file.replace(".gff", "")
                gff_path = os.path.join(gff_dir, gff_file)
                fasta_path = os.path.join(fasta_dir, f"{strain}.fasta")
                sequences = create_faidx(fasta_path, close=True)
                features = parse_gff(gff_path)
                genome_data[strain] = (sequences, features)
        
        # Get gene sequences for this cluster
        cluster_gen = iter_gene_clusters(
            genepres, genome_data, up=0, down=0, down_start_codon=False, patfilt=False, gene_list=[idx]
        )
        
        # Test cluster_cutter
        klength = 31
        stroi = {"H1-001-0005-R-G", "H1-002-0048-T-C", "H1-004-0013-A-M"}
        output_dir = os.path.join(self.test_dir, "cluster_cutter_test")
        os.makedirs(output_dir, exist_ok=True)
        
        cluster_dict_iter = []
        for cluster_gen_item in cluster_gen:
            result = cluster_cutter(
                cluster_gen_item, 
                klength, 
                stroi,
                multiple_files=False,
                canon=True,
                consider_missing_cluster=False,
                output=output_dir
            )
            cluster_dict_iter.append(result)
        
        # Check that we got results
        assert len(cluster_dict_iter) == 1
        idx_result, cluster_dict, clusterpresab_result, memchunk = cluster_dict_iter[0]
        
        assert isinstance(cluster_dict, dict)
        assert len(cluster_dict) > 0  # Should have some k-mers
        assert isinstance(memchunk, str)
        assert len(memchunk) > 0  # Should have some k-mer info

    def test_cluster_cutter_with_upstream_downstream(self):
        """Test cluster_cutter with upstream and downstream parameters"""
        # Prepare test data
        presence_absence = os.path.join(self.test_files_dir, "gene_presence_absence.csv")
        genepres = pd.read_csv(presence_absence, sep=",", index_col=0, low_memory=False)
        
        # Get first cluster
        first_cluster = genepres.iloc[0]
        idx = first_cluster.name
        clusterpresab = first_cluster.values
        
        # Prepare genome data for this cluster
        gff_dir = os.path.join(self.test_files_dir, "gffs")
        fasta_dir = os.path.join(self.test_files_dir, "fastas")
        genome_data = {}
        for gff_file in os.listdir(gff_dir):
            if gff_file.endswith(".gff"):
                strain = gff_file.replace(".gff", "")
                gff_path = os.path.join(gff_dir, gff_file)
                fasta_path = os.path.join(fasta_dir, f"{strain}.fasta")
                sequences = create_faidx(fasta_path, close=True)
                features = parse_gff(gff_path)
                genome_data[strain] = (sequences, features)
        
        # Get gene sequences for this cluster with upstream/downstream
        cluster_gen = iter_gene_clusters(
            genepres, genome_data, up=50, down=50, down_start_codon=False, patfilt=False, gene_list=[idx]
        )
        
        # Test cluster_cutter
        klength = 31
        stroi = {"H1-001-0005-R-G", "H1-002-0048-T-C", "H1-004-0013-A-M"}
        output_dir = os.path.join(self.test_dir, "cluster_cutter_updown_test")
        os.makedirs(output_dir, exist_ok=True)
        
        cluster_dict_iter = []
        for cluster_gen_item in cluster_gen:
            result = cluster_cutter(
                cluster_gen_item, 
                klength, 
                stroi,
                multiple_files=False,
                canon=True,
                consider_missing_cluster=False,
                output=output_dir
            )
            cluster_dict_iter.append(result)
        
        # Check that we got results
        assert len(cluster_dict_iter) == 1
        idx_result, cluster_dict, clusterpresab_result, memchunk = cluster_dict_iter[0]
        
        assert isinstance(cluster_dict, dict)
        assert len(cluster_dict) > 0  # Should have some k-mers
        assert isinstance(memchunk, str)
        assert len(memchunk) > 0  # Should have some k-mer info

    def test_pattern_hasher_basic(self):
        """Test pattern_hasher with basic parameters"""
        # Prepare test data
        presence_absence = os.path.join(self.test_files_dir, "gene_presence_absence.csv")
        genepres = pd.read_csv(presence_absence, sep=",", index_col=0, low_memory=False)
        
        # Prepare genome data
        gff_dir = os.path.join(self.test_files_dir, "gffs")
        fasta_dir = os.path.join(self.test_files_dir, "fastas")
        genome_data = {}
        for gff_file in os.listdir(gff_dir):
            if gff_file.endswith(".gff"):
                strain = gff_file.replace(".gff", "")
                gff_path = os.path.join(gff_dir, gff_file)
                fasta_path = os.path.join(fasta_dir, f"{strain}.fasta")
                sequences = create_faidx(fasta_path, close=True)
                features = parse_gff(gff_path)
                genome_data[strain] = (sequences, features)
        
        # Get a few clusters
        test_genes = [genepres.index[0], genepres.index[1]]
        
        # Get gene sequences for these clusters
        cluster_gen = iter_gene_clusters(
            genepres, genome_data, up=0, down=0, down_start_codon=False, patfilt=False, gene_list=test_genes
        )
        
        # Prepare cluster dict iter
        klength = 31
        stroi = {"H1-001-0005-R-G", "H1-002-0048-T-C", "H1-004-0013-A-M"}
        output_dir = os.path.join(self.test_dir, "pattern_hasher_test")
        os.makedirs(output_dir, exist_ok=True)
        
        cluster_dict_iter = []
        for cluster_gen_item in cluster_gen:
            result = cluster_cutter(
                cluster_gen_item, 
                klength, 
                stroi,
                multiple_files=False,
                canon=True,
                consider_missing_cluster=False,
                output=output_dir
            )
            cluster_dict_iter.append(result)
        
        # Test pattern_hasher
        output_dir_hasher = os.path.join(self.test_dir, "pattern_hasher_output")
        os.makedirs(output_dir_hasher, exist_ok=True)
        
        # Create output files
        hash_pat = open(os.path.join(output_dir_hasher, "hashes_to_patterns.tsv"), "w")
        kmer_hash = open(os.path.join(output_dir_hasher, "kmers_to_hashes.tsv"), "w")
        kmer_stroi = open(os.path.join(output_dir_hasher, "kmers.tsv"), "w")
        
        # Write headers
        write_headers(hash_pat, kmer_hash, genepres)
        kmer_stroi.write("cluster\tstrain\tfeature_id\tcontig\tfeature_strand\tcontig_start\tcontig_end\tgene_start\tgene_end\tstrand\tk-mer\n")
        
        # Run pattern hasher
        patterns = pattern_hasher(
            iter(cluster_dict_iter),
            kmer_stroi=kmer_stroi,
            hash_pat=hash_pat,
            kmer_hash=kmer_hash,
            genepres=genepres,
            patfilt=False,
            maf=0.05,
            output=output_dir_hasher,
            consider_missing_cluster=False
        )
        
        # Close files
        hash_pat.close()
        kmer_hash.close()
        
        # Check that patterns were generated
        assert patterns is not None
        assert len(patterns) > 0
        
        # Check that output files were created and have content
        hash_pat_file = os.path.join(output_dir_hasher, "hashes_to_patterns.tsv")
        kmer_hash_file = os.path.join(output_dir_hasher, "kmers_to_hashes.tsv")
        
        assert os.path.exists(hash_pat_file)
        assert os.path.exists(kmer_hash_file)
        
        with open(hash_pat_file, 'r') as f:
            content = f.read()
            assert len(content) > 0
            # Should have header + at least one pattern
            lines = content.strip().split('\n')
            assert len(lines) >= 2
        
        with open(kmer_hash_file, 'r') as f:
            content = f.read()
            assert len(content) > 0
            # Should have header + at least one kmer
            lines = content.strip().split('\n')
            assert len(lines) >= 2

    def test_pattern_hasher_with_maf_filtering(self):
        """Test pattern_hasher with MAF filtering"""
        # Prepare test data
        presence_absence = os.path.join(self.test_files_dir, "gene_presence_absence.csv")
        genepres = pd.read_csv(presence_absence, sep=",", index_col=0, low_memory=False)
        
        # Prepare genome data
        gff_dir = os.path.join(self.test_files_dir, "gffs")
        fasta_dir = os.path.join(self.test_files_dir, "fastas")
        genome_data = {}
        for gff_file in os.listdir(gff_dir):
            if gff_file.endswith(".gff"):
                strain = gff_file.replace(".gff", "")
                gff_path = os.path.join(gff_dir, gff_file)
                fasta_path = os.path.join(fasta_dir, f"{strain}.fasta")
                sequences = create_faidx(fasta_path, close=True)
                features = parse_gff(gff_path)
                genome_data[strain] = (sequences, features)
        
        # Get a cluster
        test_genes = [genepres.index[0]]
        
        # Get gene sequences for this cluster
        cluster_gen = iter_gene_clusters(
            genepres, genome_data, up=0, down=0, down_start_codon=False, patfilt=False, gene_list=test_genes
        )
        
        # Prepare cluster dict iter
        klength = 31
        stroi = {"H1-001-0005-R-G", "H1-002-0048-T-C", "H1-004-0013-A-M"}
        output_dir = os.path.join(self.test_dir, "pattern_hasher_maf_test")
        os.makedirs(output_dir, exist_ok=True)
        
        cluster_dict_iter = []
        for cluster_gen_item in cluster_gen:
            result = cluster_cutter(
                cluster_gen_item, 
                klength, 
                stroi,
                multiple_files=False,
                canon=True,
                consider_missing_cluster=False,
                output=output_dir
            )
            cluster_dict_iter.append(result)
        
        # Test pattern_hasher with high MAF threshold
        output_dir_hasher = os.path.join(self.test_dir, "pattern_hasher_maf_output")
        os.makedirs(output_dir_hasher, exist_ok=True)
        
        # Create output files
        hash_pat = open(os.path.join(output_dir_hasher, "hashes_to_patterns.tsv"), "w")
        kmer_hash = open(os.path.join(output_dir_hasher, "kmers_to_hashes.tsv"), "w")
        kmer_stroi = open(os.path.join(output_dir_hasher, "kmers.tsv"), "w")
        
        # Write headers
        write_headers(hash_pat, kmer_hash, genepres)
        kmer_stroi.write("cluster\tstrain\tfeature_id\tcontig\tfeature_strand\tcontig_start\tcontig_end\tgene_start\tgene_end\tstrand\tk-mer\n")
        
        # Run pattern hasher with high MAF (should filter out more patterns)
        patterns_high_maf = pattern_hasher(
            iter(cluster_dict_iter),
            kmer_stroi=kmer_stroi,
            hash_pat=hash_pat,
            kmer_hash=kmer_hash,
            genepres=genepres,
            patfilt=False,
            maf=0.5,  # High MAF threshold
            output=output_dir_hasher,
            consider_missing_cluster=False
        )
        
        # Close files
        hash_pat.close()
        kmer_hash.close()
        
        # Check that patterns were generated (but likely fewer than with lower MAF)
        assert patterns_high_maf is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
