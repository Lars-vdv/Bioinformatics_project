from bio_seq import bio_seq

revp = bio_seq()
revp = parse_fasta("test_data/rosalind_revp.txt", seq_type="DNA", create_variables=True)

revp_rc = getattr(seq_1, "seq")