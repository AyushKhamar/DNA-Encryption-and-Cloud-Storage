# encryption code
import string
from time import time
import utils
from utils import *

# variables that are used everywhere
number_of_rounds = None
chromosome_length = None
decryption_key = None

# use all the functions that are used for identification in decryption
identifier_key = key_del
identifier_no_rounds = no_rounds_del
identifier_round = round_del
identifier_reshape = reshape_del
identifier_crossover = crossover_del
identifier_crossover_type = crossover_type_del
identifier_single_point_crossover = single_point_crossover_del
identifier_rotate_crossover = rotate_crossover_del
identifier_rotation_offset = rotation_offset_del
identifier_rotation_types = rotation_types_del
identifier_mutation = mutation_del
identifier_complement_mutation = complement_mutation_del
identifier_alter_mutation = alter_mutation_del
identifier_mutation_table = mutation_table_del
identifier_chromosome = chromosome_del


def initialize_globals():
    global number_of_rounds
    global decryption_key
    number_of_rounds = random.randrange(3, 12, 2)
    decryption_key = ""


def encrypt_key(plaintext, key):
    # we will xor the data with to key for first level encruption
    if len(plaintext) > len(key):
        factor = int(len(plaintext) / len(key))
        key += key * factor
        return bitxor(plaintext, key)
    return bitxor(plaintext, key)


def reshape_of_dna(dna_sequence):
    # here we will make the population through a single individual, we will make an array of the population generated and return the same
    global chromosome_length
    global decryption_key
    # this is for generation of the population
    divs = divisors(len(dna_sequence))
    total_chromosomes = divs[random.randint(0, len(divs) - 1)]
    chromosome_length = int(len(dna_sequence) / total_chromosomes)
    chromosomes = []
    decryption_key += identifier_reshape + \
        str(chromosome_length) + identifier_reshape
    # this is taking out the created population
    for i in range(0, len(dna_sequence), chromosome_length):
        chromosomes.append(dna_sequence[i:i + chromosome_length])

    return chromosomes


def reverse_reshape(total_population):
    # convert the chromosome total_population back to DNA sequence
    return "".join(total_population)


def crossover_using_rotate(total_population):
    # we have introduced random probability and based on the results we will choose whether to rotate the current chromosome left or right
    global chromosome_length
    global decryption_key
    newly_generated_total_population = []
    decryption_key += identifier_rotate_crossover
    # we have already defined a rotation value which we will vary with every round of encryption
    rotation_offset = random.randint(1, chromosome_length)
    decryption_key += identifier_rotation_offset + \
        str(rotation_offset) + identifier_rotation_offset
    decryption_key += identifier_rotation_types

    for chromosome in total_population:
        probability = random.uniform(0, 1)
        if probability > 0.5:
            decryption_key += "right|"
            right_part_one = chromosome[0: len(chromosome) - rotation_offset]
            right_part_two = chromosome[len(chromosome) - rotation_offset:]
            newly_generated_total_population.append(
                right_part_two + right_part_one)
        else:
            decryption_key += "left|"
            left_part_one = chromosome[0: rotation_offset]
            left_part_two = chromosome[rotation_offset:]
            newly_generated_total_population.append(
                left_part_two + left_part_one)

    decryption_key += identifier_rotation_types
    decryption_key += identifier_rotate_crossover

    return newly_generated_total_population


def crossover_using_single_point(total_population):
    # we will use the single point crossover method to choose two parents and make their offspring
    global decryption_key
    decryption_key += identifier_single_point_crossover

    newly_generated_total_population = []
    for i in range(0, len(total_population) - 1, 2):
        candidate1 = total_population[i]
        candidate2 = total_population[i + 1]
        # we know that chromosomes will have the same length and now we generate the point of crossover using random function
        length = len(candidate1)
        crossover_point = random.randint(0, length - 1)
        decryption_key += str(crossover_point) + "|"
        first_offspring = candidate2[0: crossover_point] + \
            candidate1[crossover_point:]
        second_offspring = candidate1[0: crossover_point] + \
            candidate2[crossover_point:]
        newly_generated_total_population.append(first_offspring)
        newly_generated_total_population.append(second_offspring)

    if len(total_population) % 2 == 1:
        newly_generated_total_population.append(
            total_population[len(total_population) - 1])

    decryption_key += identifier_single_point_crossover
    return newly_generated_total_population


def crossover(total_population):
    global decryption_key
    # now we will choose what kind of crossover to apply to a certain population based on probability
    probability = random.uniform(0, 1)

    if probability < 0.33:
        decryption_key += identifier_crossover_type + \
            "rotate_crossover" + identifier_crossover_type
        return crossover_using_rotate(total_population)
    elif probability >= 0.33 and probability < 0.66:
        decryption_key += identifier_crossover_type + \
            "single_point_crossover" + identifier_crossover_type
        return crossover_using_single_point(total_population)
    else:
        decryption_key += identifier_crossover_type + "both" + identifier_crossover_type
        total_population = crossover_using_rotate(total_population)
        return crossover_using_single_point(total_population)


def chromosome_complement(chromosome, p1, p2):
    # the basic idea is to invert all bits between the two points p1 and p2
    newly_generated_chromosome = ""
    for i in range(len(chromosome)):
        if i >= p1 and i <= p2:
            if chromosome[i] == '0':
                newly_generated_chromosome += '1'
            else:
                newly_generated_chromosome += '0'
        else:
            newly_generated_chromosome += chromosome[i]

    return newly_generated_chromosome


def alter_dna_bases(current_bases):
    # we know that A pairs with T and C pairs with G, so we will change the DNA bases from C --> G and A --> T randomly
    alter_dna_table = {}

    for _ in range(2):
        # it basically works like a swap operation, we take one and then remove it from the list for the two bases wehave and then we swap those components
        base1 = current_bases[random.randint(0, len(current_bases) - 1)]
        current_bases.remove(base1)
        # choose one randomly then remove it from list
        base2 = current_bases[random.randint(0, len(current_bases) - 1)]
        current_bases.remove(base2)
        # now we swap the two
        alter_dna_table[base1] = base2
        alter_dna_table[base2] = base1

    return alter_dna_table


def mutation(total_population):
    # mutation is to ensure that we dont get the same offspring again and again, it is forced evolution, we will combine the complement and change dna functions in order to give out a mutation
    global decryption_key

    bases = ['A', 'C', 'G', 'T']
    alter_dna_table = alter_dna_bases(bases)
    decryption_key += identifier_mutation_table + \
        str(alter_dna_table) + identifier_mutation_table

    newly_generated_total_population = []
    for chromosome in total_population:
        decryption_key += identifier_chromosome
        # first operate on complement
        temp_chromosome = dna_to_bits(
            chromosome, utils.dna_base_to_two_bits_table)
        decryption_key += identifier_complement_mutation
        p1 = random.randint(0, len(temp_chromosome) - 1)
        p2 = random.randint(p1, len(temp_chromosome) - 1)
        decryption_key += "(%s, %s)" % (p1, p2)
        decryption_key += identifier_complement_mutation
        temp_chromosome = chromosome_complement(temp_chromosome, p1, p2)
        # well convert to dna base table
        temp_vector = group_bits(temp_chromosome, 4)
        last_dna_base = None
        # now we will convert only if the last element is not of length 2
        if len(temp_vector[len(temp_vector) - 1]) == 2:
            last_dna_base = utils.two_bits_to_dna_base_table[temp_vector[len(
                temp_vector) - 1]]
            # dont convert the 2 bit elements
            temp_vector = temp_vector[:-1]

        newly_generated_dna_sequence = bits_to_dna(
            temp_vector, utils.four_bits_to_two_dna_base_table)
        if last_dna_base is not None:
            newly_generated_dna_sequence += last_dna_base

        # and then alter the dna bases between p1 and p2
        decryption_key += identifier_alter_mutation
        p1 = random.randint(0, len(newly_generated_dna_sequence) - 1)
        p2 = random.randint(p1, len(newly_generated_dna_sequence) - 1)
        decryption_key += "(%s, %s)" % (p1, p2)
        decryption_key += identifier_alter_mutation
        newly_generated_chromosome = ""
        for i in range(len(newly_generated_dna_sequence)):
            if i >= p1 and i <= p2:
                newly_generated_chromosome += alter_dna_table[newly_generated_dna_sequence[i]]
            else:
                newly_generated_chromosome += newly_generated_dna_sequence[i]

        newly_generated_total_population.append(newly_generated_chromosome)
        decryption_key += identifier_chromosome

    return newly_generated_total_population


def dna_get(text, key):
    global number_of_rounds
    global decryption_key
    # we first need to convert the data to binary to work on it
    binary_data_1 = binarized_data(text)
    dna_seq = bits_to_dna(binary_data_1, utils.two_bits_to_dna_base_table)
    binary_data_2 = dna_seq
    decryption_key += identifier_no_rounds + \
        str(number_of_rounds) + identifier_no_rounds
    # now we will iterate number_of_rounds time throught the loop to rehape crossover and mutate
    while number_of_rounds > 0:
        decryption_key += identifier_round
        # we will encrypt the data with the key we have and then we can convert the binary data back to a dna sequence
        binary_data_2 = bits_to_dna(
            group_bits(encrypt_key(dna_to_bits(reverse_reshape(
                binary_data_2), utils.dna_base_to_two_bits_table), key)),
            utils.two_bits_to_dna_base_table)
        # we can move ahead to create the individual population now
        binary_data_2 = reshape_of_dna(binary_data_2)
        # now we will pick individuals and use probability to apply crossover to them        decryption_key += identifier_crossover
        binary_data_2 = crossover(binary_data_2)
        decryption_key += identifier_crossover
        # now to remove the possibility of all offsprings becoming the same we apply mutation to keep the evolution going on
        decryption_key += identifier_mutation
        binary_data_2 = mutation(binary_data_2)
        decryption_key += identifier_mutation

        number_of_rounds -= 1
        decryption_key += identifier_round

    return reverse_reshape(binary_data_2)


def main(original_filename):
    global decryption_key

    original_file = open(original_filename, "r")
    plaintext = original_file.read()
    # we generate a random 128 bit key length key for encryption
    key = str2bin(''.join(random.SystemRandom().choice(
        string.ascii_letters + string.digits) for _ in range(16)))

    initialize_globals()
    decryption_key += identifier_key + key + identifier_key
    # we need the encoding tables to move ahead
    generate_pre_processing_tables()
    generate_mutation_tables()
    # now we can encrypt the plaintext that we extracted from the plain file
    encrypted_text = dna_get(plaintext, key)
    key_file = open(key_filename, "w")
    encrypted_file = open(encrypted_filename, "w")
    encrypted_file.write(encrypted_text)
    key_file.write(decryption_key)

    encrypted_file.close()
    original_file.close()
    key_file.close()
    return encrypted_text
