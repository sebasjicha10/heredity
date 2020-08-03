import csv
import itertools
import sys
import random
import math


PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """

    joint_probabilities = []

    # Define the gene probabilities for people with no parents
    def no_parent_probability(person, genes):
        return PROBS["gene"][genes]

    # Define the gene probability for people with parents 
    def gene_probability_parents(person, genes):

        # Variables
        mother = people[person]["mother"]
        father = people[person]["father"]
        mother_probability = 0
        father_probability = 0

        # Probabilities for each parent
        if mother in one_gene:
            mother_probability = 0.5 
   
        elif mother in two_genes:
            mother_probability = 1 - PROBS["mutation"]

        else:
            mother_probability = PROBS["mutation"]

        if father in one_gene:
            father_probability = 0.5 
            
        elif father in two_genes:
            father_probability = 1 - PROBS["mutation"]

        else:
            father_probability = PROBS["mutation"]


        # Probability of 0 copies of the gene
        if genes == 0:
            return ((1 - mother_probability) * (1 - father_probability))

        # Probability of 1 copy of the gene
        elif genes == 1:
            return ((1 - mother_probability) * father_probability) + (mother_probability * (1 - father_probability))

        # Probability of 2 copies of the gene
        elif genes == 2:
            return (mother_probability * father_probability)

   
    def get_trait_probability(person, genes):

        # People with trait
        if person in have_trait:
            return PROBS["trait"][genes][True]

        else:
            return PROBS["trait"][genes][False]

    # Iterate over all people
    for person in people:

        gene_probability = 0
        trait_probability = 0
        
        # Check how many genes to test for
        if person in one_gene:
            if not people[person]["mother"]:
                gene_probability = no_parent_probability(person, 1)

            else:
                gene_probability = gene_probability_parents(person, 1)
            trait_probability = get_trait_probability(person, 1)

        elif person in two_genes:
            if not people[person]["mother"]:
                gene_probability = no_parent_probability(person, 2)

            else:
                gene_probability = gene_probability_parents(person, 2)
            trait_probability = get_trait_probability(person, 2)

        # Probability of person having 0 genes
        else:
            if not people[person]["mother"]:
                gene_probability = no_parent_probability(person, 0)

            else:
                gene_probability = gene_probability_parents(person, 0)
            trait_probability = get_trait_probability(person, 0)
            
        final_individual_probability = gene_probability * trait_probability

        joint_probabilities.append(final_individual_probability)

    return math.prod(joint_probabilities)


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """

    for person in probabilities:

        # Add gene probility
        genes = 0

        if person in one_gene:
            genes = 1

        elif person in two_genes:
            genes = 2

        probabilities[person]["gene"][genes] = probabilities[person]["gene"][genes] + p

        # Add trait probability
        trait = False

        if person in have_trait:
            trait = True

        probabilities[person]["trait"][trait] = probabilities[person]["trait"][trait] + p


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """

    for person in probabilities:

        # Normalize gene
        gene_2 = probabilities[person]["gene"][2]
        gene_1 = probabilities[person]["gene"][1]
        gene_0 = probabilities[person]["gene"][0]

        gene_normalizer = 1 / (gene_0 + gene_1 + gene_2)

        # Update genes
        probabilities[person]["gene"][2] = gene_2 * gene_normalizer
        probabilities[person]["gene"][1] = gene_1 * gene_normalizer
        probabilities[person]["gene"][0] = gene_0 * gene_normalizer

        # Normalize trait
        true_trait = probabilities[person]["trait"][True]
        false_trait = probabilities[person]["trait"][False]

        trait_normalizer = 1 / (true_trait + false_trait)

        # Update traits
        probabilities[person]["trait"][True] = true_trait * trait_normalizer
        probabilities[person]["trait"][False] = false_trait * trait_normalizer


if __name__ == "__main__":
    main()
