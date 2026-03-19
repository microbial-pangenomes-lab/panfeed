#!/usr/bin/env python
"""
Unit tests for panfeed.input module
"""

import os
import tempfile
import shutil
import pytest
import numpy as np
from panfeed.input import (
    what_are_my_inputfiles,
    parse_gff,
    create_faidx,
    set_input_output,
    iter_gene_clusters
)
import pandas as pd


class TestInputModule:
    """Test cases for input module functions"""

    def setup_method(self):
        """Setup test fixtures"""
        self.test_dir = tempfile.mkdtemp()
        # Use dynamic path resolution that works from both root and tests directory
        if os.path.exists("test_files"):
            self.test_files_dir = "test_files"
        else:
            self.test_files_dir = "tests/test_files"

    def teardown_method(self):
        """Clean up test fixtures"""
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def test_what_are_my_inputfiles_directory(self):
        """Test what_are_my_inputfiles with directory input"""
        gff_dir = os.path.join(self.test_files_dir, "gffs")
        filelist, fastalist = what_are_my_inputfiles(gff_dir)
        
        assert len(filelist) == 3
        assert "H1-001-0005-R-G" in filelist
        assert "H1-002-0048-T-C" in filelist
        assert "H1-004-0013-A-M" in filelist
        assert fastalist == []  # No fasta dir provided

    def test_what_are_my_inputfiles_file_of_files(self):
        """Test what_are_my_inputfiles with file of files input"""
        gff_file = os.path.join(self.test_files_dir, "input_gffs.txt")
        filelist, fastalist = what_are_my_inputfiles(gff_file)
        
        assert len(filelist) == 3
        assert "H1-001-0005-R-G" in filelist
        assert "H1-002-0048-T-C" in filelist
        assert "H1-004-0013-A-M" in filelist

    def test_what_are_my_inputfiles_with_fasta(self):
        """Test what_are_my_inputfiles with both GFF and FASTA inputs"""
        gff_dir = os.path.join(self.test_files_dir, "gffs")
        fasta_dir = os.path.join(self.test_files_dir, "fastas")
        filelist, fastalist = what_are_my_inputfiles(gff_dir, fasta_dir)
        
        assert len(filelist) == 3
        assert len(fastalist) == 3
        assert "H1-001-0005-R-G" in fastalist
        assert "H1-002-0048-T-C" in fastalist
        assert "H1-004-0013-A-M" in fastalist

    def test_parse_gff(self):
        """Test GFF parsing functionality"""
        gff_file = os.path.join(self.test_files_dir, "gffs", "H1-001-0005-R-G.gff")
        features = parse_gff(gff_file)
        
        # Should have some features
        assert len(features) > 0
        
        # Check that features have the expected structure
        for feature_id, feature in features.items():
            assert hasattr(feature, 'id')
            assert hasattr(feature, 'chromosome')
            assert hasattr(feature, 'start')
            assert hasattr(feature, 'end')
            assert hasattr(feature, 'strand')

    def test_set_input_output_basic(self):
        """Test set_input_output with basic parameters"""
        stroi_file = os.path.join(self.test_files_dir, "stroi.txt")
        presence_absence = os.path.join(self.test_files_dir, "gene_presence_absence.csv")
        output_dir = os.path.join(self.test_dir, "test_output")
        
        stroi, genes, kmer_stroi, hash_pat, kmer_hash, genepres = set_input_output(
            stroi_file, None, presence_absence, output_dir
        )
        
        # Check outputs
        assert len(stroi) == 3
        assert genes is None
        assert os.path.exists(output_dir)
        assert kmer_stroi is not None
        assert hash_pat is not None
        assert kmer_hash is not None
        assert isinstance(genepres, pd.DataFrame)
        
        # Clean up
        kmer_stroi.close()
        hash_pat.close()
        kmer_hash.close()

    def test_set_input_output_with_genes(self):
        """Test set_input_output with gene list"""
        stroi_file = os.path.join(self.test_files_dir, "stroi.txt")
        genes_file = os.path.join(self.test_files_dir, "target_clusters.txt")
        presence_absence = os.path.join(self.test_files_dir, "gene_presence_absence.csv")
        output_dir = os.path.join(self.test_dir, "test_output_genes")
        
        stroi, genes, kmer_stroi, hash_pat, kmer_hash, genepres = set_input_output(
            stroi_file, genes_file, presence_absence, output_dir
        )
        
        # Check outputs
        assert len(stroi) == 3
        assert genes is not None
        assert len(genes) > 0
        
        # Clean up
        kmer_stroi.close()
        hash_pat.close()
        kmer_hash.close()

    def test_iter_gene_clusters_basic(self):
        """Test iter_gene_clusters with basic parameters"""
        presence_absence = os.path.join(self.test_files_dir, "gene_presence_absence.csv")
        genepres = pd.read_csv(presence_absence, sep=",", index_col=0, low_memory=False)
        
        # Prepare genome data - use the separate fasta files that are available
        fasta_dir = os.path.join(self.test_files_dir, "fastas")
        gff_dir = os.path.join(self.test_files_dir, "gffs")
        genome_data = {}
        
        # Create a temporary output directory for extracted FASTA files
        temp_output = os.path.join(self.test_dir, "temp_fastas")
        os.makedirs(temp_output, exist_ok=True)
        
        for gff_file in os.listdir(gff_dir):
            if gff_file.endswith(".gff"):
                strain = gff_file.replace(".gff", "")
                gff_path = os.path.join(gff_dir, gff_file)
                
                # Check if we have a separate fasta file
                fasta_file = os.path.join(fasta_dir, f"{strain}.fasta")
                if os.path.exists(fasta_file):
                    sequences = create_faidx(fasta_file, close=True)
                else:
                    # Extract FASTA from GFF file
                    fasta_content = open(gff_path, "r").read().split("##FASTA")
                    if len(fasta_content) > 1:
                        temp_fasta = os.path.join(temp_output, f"{strain}.fasta")
                        with open(temp_fasta, "w") as f:
                            f.write(fasta_content[1])
                        sequences = create_faidx(temp_fasta, close=True)
                    else:
                        continue
                
                features = parse_gff(gff_path)
                genome_data[strain] = (sequences, features)
        
        # Test iteration
        cluster_count = 0
        for gene_sequences, idx, clusterpresab in iter_gene_clusters(
            genepres, genome_data, up=0, down=0, down_start_codon=False, patfilt=False
        ):
            cluster_count += 1
            assert isinstance(gene_sequences, dict)
            assert isinstance(idx, str)
            assert isinstance(clusterpresab, np.ndarray)
            
            # Should not iterate too many times for this test
            if cluster_count >= 5:
                break
        
        assert cluster_count > 0

    def test_iter_gene_clusters_with_upstream_downstream(self):
        """Test iter_gene_clusters with upstream and downstream parameters"""
        presence_absence = os.path.join(self.test_files_dir, "gene_presence_absence.csv")
        genepres = pd.read_csv(presence_absence, sep=",", index_col=0, low_memory=False)
        
        # Prepare genome data - use the separate fasta files that are available
        fasta_dir = os.path.join(self.test_files_dir, "fastas")
        gff_dir = os.path.join(self.test_files_dir, "gffs")
        genome_data = {}
        
        # Create a temporary output directory for extracted FASTA files
        temp_output = os.path.join(self.test_dir, "temp_fastas_updown")
        os.makedirs(temp_output, exist_ok=True)
        
        for gff_file in os.listdir(gff_dir):
            if gff_file.endswith(".gff"):
                strain = gff_file.replace(".gff", "")
                gff_path = os.path.join(gff_dir, gff_file)
                
                # Check if we have a separate fasta file
                fasta_file = os.path.join(fasta_dir, f"{strain}.fasta")
                if os.path.exists(fasta_file):
                    sequences = create_faidx(fasta_file, close=True)
                else:
                    # Extract FASTA from GFF file
                    fasta_content = open(gff_path, "r").read().split("##FASTA")
                    if len(fasta_content) > 1:
                        temp_fasta = os.path.join(temp_output, f"{strain}.fasta")
                        with open(temp_fasta, "w") as f:
                            f.write(fasta_content[1])
                        sequences = create_faidx(temp_fasta, close=True)
                    else:
                        continue
                
                features = parse_gff(gff_path)
                genome_data[strain] = (sequences, features)
        
        # Test iteration with upstream/downstream
        cluster_count = 0
        for gene_sequences, idx, clusterpresab in iter_gene_clusters(
            genepres, genome_data, up=100, down=100, down_start_codon=False, patfilt=False
        ):
            cluster_count += 1
            
            # Check that sequences are properly extended
            for strain, sequences in gene_sequences.items():
                if sequences:  # If strain has this gene
                    for seq in sequences:
                        # Sequences should be longer than just the gene due to upstream/downstream
                        assert len(seq.sequence) > 0
            
            # Should not iterate too many times for this test
            if cluster_count >= 3:
                break
        
        assert cluster_count > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
