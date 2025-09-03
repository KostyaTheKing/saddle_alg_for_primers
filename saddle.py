import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import copy


class SaddleSimulator:
    """
    A class that chooses a set of primers with the SADDLE algorithm.
    """

    def __init__(self, all_primers: pd.DataFrame, gt=100, len_substring=4, seed=None):
        self.all_primers = all_primers
        self.gt = gt
        self.len_substring = len_substring
        self.rng = np.random.default_rng(seed)
        self.loss_progress = []
        self.primer_variants = None

    def _generate_random_set(self, current_primer_list=[]) -> list[str]:
        """
        Generates a random list of primers or replaces a random primer with a new one.
        """
        # Get number of columns
        _, m = self.all_primers.shape

        if not current_primer_list:
            # Initialize empty list of random primers
            random_primer_list = []

            # Add random primer from each column to the list
            for i in range(m):
                random_primer = (
                    self.all_primers.iloc[:, i]
                    .dropna()
                    .sample(n=1, random_state=self.rng)
                    .to_list()
                )
                random_primer_list += random_primer
            return random_primer_list
        else:
            # Check whether a primer can be changed
            new_primer_list = copy.deepcopy(current_primer_list)
            series_len = 0
            while series_len == 0:
                # Choose a primer to change in a current list
                index_of_a_primer_to_change = self.rng.integers(0, m)
                current_primer = new_primer_list[index_of_a_primer_to_change]
                # Get a new non duplicate primer
                primer_series = (
                    self.all_primers.iloc[:, index_of_a_primer_to_change]
                    .dropna()
                    .reset_index(drop=True)
                )
                series_len = len(primer_series[primer_series != current_primer])
            # Get a new primer
            new_primer = (
                primer_series[primer_series != current_primer]
                .sample(1, random_state=self.rng)
                .to_list()[0]
            )
            # Replace an old primer with a new primer
            new_primer_list[index_of_a_primer_to_change] = new_primer
            return new_primer_list

    @staticmethod
    def _distance(primer: str, substring: str) -> int:
        """
        Calculates a distance from 3'-end of a primer to a substring.
        """
        substring = substring[::-1]
        primer = primer[::-1]

        return primer.index(substring)

    @staticmethod
    def _revcomp(dna_seq: str) -> str:
        """
        Returns a reverse-complement sequence of a given DNA sequence.
        """
        # Watson–Crick base pairs
        nucl_complement = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
        }

        return "".join(nucl_complement[i] for i in list(dna_seq[::-1]))

    def _loss_calc(self, current_primer_list: list) -> float:
        """
        A loss function.
        """
        # Creating an empty list.
        H = list()
        # Iterate over all primers in a list
        for i in current_primer_list:
            # Get last len_substring nucleotides of a primer
            i = i[-self.len_substring :]
            i = i.upper()
            # Count GC content
            numGC = i.count("G") + i.count("C")
            # Count how many times i is found in a primer list
            for_numb = sum(1 for j in current_primer_list if i in j)
            # Calculate a score for a loss function
            score_1 = sum(
                1 / (self._distance(j, i) + 1) for j in current_primer_list if i in j
            )
            # Get a reverse-complement of the substring
            rev_comp_i = self._revcomp(i)
            # Count how many times rev_comp_i is found in a primer list
            rev_numb = sum(1 for j in current_primer_list if rev_comp_i in j)
            # Calculate a score for a loss function
            score_2 = sum(
                1 / (self._distance(j, rev_comp_i) + 1)
                for j in current_primer_list
                if rev_comp_i in j
            )
            # Add calculated variables to the list
            H.append((numGC, for_numb, score_1, rev_numb, score_2))
        L = sum(
            2 ** (self.len_substring + i[0]) * (i[1] * i[3]) * i[2] * i[4] for i in H
        )
        return L

    def saddle(self):
        """
        SADDLE algorithm. This function returns:
        1) The last list of primers that it generated.
        2) Adds a loss progress to a class instance.
        3) Adds a list of primers generated at each stage to a class instance.
        """
        # Total amount of generation. After gt SADDLE is changed to SGD.
        last_generation = int(self.gt * 1.5)

        # List of calculated loss at each stage
        loss_progress = []

        # All primer variants selected at each stage
        primer_variants = []

        # Creating initial primer set
        sg = self._generate_random_set()

        primer_variants.append(sg)

        for generation in range(last_generation):
            # Loss of the old primer set
            L_Sg = self._loss_calc(sg)
            # New primer set
            t = self._generate_random_set(sg)
            # Loss of a new primer set
            L_T = self._loss_calc(t)

            if L_T < L_Sg:
                sg = copy.deepcopy(t)
            else:
                if generation < self.gt:
                    prob = np.exp((L_Sg - L_T) / (1 / (generation + 1)))
                    if self.rng.uniform(0, 1) < prob:
                        sg = copy.deepcopy(t)
            # Adding loss information to a list
            loss_progress.append(L_Sg)
            primer_variants.append(sg)

        self.loss_progress = np.array([np.arange(0, last_generation, 1), loss_progress])
        self.primer_variants = pd.DataFrame(
            primer_variants, columns=self.all_primers.columns
        )
        return sg

    def plot_loss(self, **kwargs) -> None:
        '''
        Make a plot of loss to generation stage.
        '''
        plt.plot(self.loss_progress[0], self.loss_progress[1], **kwargs)