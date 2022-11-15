## Task 1.a

Download muscle from [here](https://drive5.com/muscle5/).

Create multiple sequence alignment:

```{bash}
path_to_Muscle -super5 human.fa -output human_aligned.fasta
path_to_Muscle -super5 mouse.fa -output mouse_aligned.fasta
```

## Task1.b

Install weblogo with `pip` or `conda`. Then:

```{bash}
weblogo --format PNG --size large < human_aligned.fasta > h_logo.png
weblogo --format PNG --size large < mouse_aligned.fasta > m_logo.png
```

