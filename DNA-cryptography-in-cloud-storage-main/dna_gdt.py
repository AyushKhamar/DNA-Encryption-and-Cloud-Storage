# decryption code
from time import time
import ast
import utils
from utils import *

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

def encryption_of_key(plaintext, key):
    # we will xor the data with to key for first level encruption
    if len(plaintext) > len(key):
        factor = int(len(plaintext) / len(key))
        key += key * factor
        return bitxor(plaintext, key)
    return bitxor(plaintext, key)


def reshape_of_dna(dna_sequence, reshape_info):
    # here we will make the population through a single individual, we will make an array of the population generated and return the same
    chromosome_length = int(reshape_info[0])
    newly_generated_chromosomes = []
    for i in range(0, len(dna_sequence), chromosome_length):
        newly_generated_chromosomes.append(
            dna_sequence[i:i + chromosome_length])

    return newly_generated_chromosomes


def reverse_reshape(total_population):
    # convert the chromosome population back to DNA sequence
    return "".join(total_population)


def crossover_using_rotate(total_population, rotate_info):
    # we have introduced random probability and based on the results we will choose whether to rotate the current chromosome left or right
    newly_generated_population = []
    # we have already defined a rotation value which we will vary with every round of encryption
    rotation_offset = int(get_pattern(identifier_rotation_offset, rotate_info)[0])
    rotations = get_pattern(identifier_rotation_types, rotate_info)[0].split("|")[:-1]

    for i in range(len(total_population)):
        chromosome = total_population[i]
        direction = rotations[i]

        if direction == "left":
            right_first = chromosome[0: len(chromosome) - rotation_offset]
            right_second = chromosome[len(chromosome) - rotation_offset:]
            newly_generated_population.append(right_second + right_first)
        elif direction == "right":
            left_first = chromosome[0: rotation_offset]
            left_second = chromosome[rotation_offset:]
            newly_generated_population.append(left_second + left_first)

    return newly_generated_population


def crossover_using_single_point(total_population, single_point_info):
    # we will use the single point crossover method to choose two parents and make their offspring
    crossover_points = [int(p)
                        for p in single_point_info.split("|") if p != '']

    newly_generated_population = []
    for i in range(0, len(total_population) - 1, 2):
        parent_one = total_population[i]
        parent_two = total_population[i + 1]
        # we know that chromosomes will have the same length and now we generate the point of crossover using random function
        crossover_point = crossover_points[int(i / 2)]
        first_offspring = parent_two[0: crossover_point] + \
            parent_one[crossover_point:]
        second_offspring = parent_one[0: crossover_point] + \
            parent_two[crossover_point:]
        newly_generated_population.append(first_offspring)
        newly_generated_population.append(second_offspring)

    if len(total_population) % 2 == 1:
        newly_generated_population.append(
            total_population[len(total_population) - 1])

    return newly_generated_population


def crossover(total_population, crossover_info):
    # now we will choose what kind of crossover to apply to a certain population based on probability
    crossover_type = get_pattern(identifier_crossover_type, crossover_info)[0]

    if crossover_type == "rotate_crossover":
        rotate_info = get_pattern(identifier_rotate_crossover, crossover_info)[0]
        return crossover_using_rotate(total_population, rotate_info)
    elif crossover_type == "single_point_crossover":
        single_point_info = get_pattern(
            identifier_single_point_crossover, crossover_info)[0]
        return crossover_using_single_point(total_population, single_point_info)
    elif crossover_type == "both":
        rotate_info = get_pattern(identifier_rotate_crossover, crossover_info)[0]
        single_point_info = get_pattern(
            identifier_single_point_crossover, crossover_info)[0]
        total_population = crossover_using_single_point(
            total_population, single_point_info)
        return crossover_using_rotate(total_population, rotate_info)


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


def mutation(total_population, mutation_info):
    # mutation is to ensure that we dont get the same offspring again and again, it is forced evolution, we will combine the complement and change dna functions in order to give out a mutation
    alter_dna_table = ast.literal_eval(
        get_pattern(identifier_mutation_table, mutation_info[0])[0])

    chromosomes_info = get_pattern(identifier_chromosome, mutation_info[0])

    newly_generated_population = []
    for i in range(len(total_population)):
        chromosome = total_population[i]
        chromosome_info = chromosomes_info[i]
        # we use points p1 and p2 to alter the dna bases
        alter_info = get_pattern(identifier_alter_mutation, chromosome_info)[0]
        p1, p2 = ast.literal_eval(alter_info)
        newly_generated_chromosome = ""
        for i in range(len(chromosome)):
            if i >= p1 and i <= p2:
                newly_generated_chromosome += alter_dna_table[chromosome[i]]
            else:
                newly_generated_chromosome += chromosome[i]

        temp_vector = group_bases(newly_generated_chromosome)
        last_two_bits = None
        if len(newly_generated_chromosome) % 2 == 1:
            last_two_bits = utils.dna_base_to_two_bits_table[newly_generated_chromosome[-1]]

            temp_vector = temp_vector[:-1]

        bits_seq = dna_to_bits(
            temp_vector, utils.two_dna_base_to_four_bits_table)

        if last_two_bits is not None:
            bits_seq += last_two_bits

        complement_info = get_pattern(
            identifier_complement_mutation, chromosome_info)[0]
        p1, p2 = ast.literal_eval(complement_info)
        temp_chromosome = chromosome_complement(bits_seq, p1, p2)
        temp_chromosome = group_bits(temp_chromosome)
        newly_generated_chromosome = bits_to_dna(
            temp_chromosome, utils.two_bits_to_dna_base_table)

        newly_generated_population.append(newly_generated_chromosome)

    return newly_generated_population


def dna_gdt(plaintext, key):
    # we first need to convert the data to binary to work on it
    number_of_rounds = int(get_pattern(identifier_no_rounds, key)[0])
    rounds = get_pattern(identifier_round, key)

    binary_sequence = plaintext
    # now we will iterate number_of_rounds time throught the loop to rehape crossover and mutate
    # run the algorithm "number_of_rounds" times
    while number_of_rounds > 0:
        round_info = rounds[number_of_rounds - 1]

        # generation of individual population
        binary_sequence = reshape_of_dna(
            binary_sequence, get_pattern(identifier_reshape, round_info))

        # mutation
        binary_sequence = mutation(
            binary_sequence, get_pattern(identifier_mutation, round_info))

        # crossover based on probability
        binary_sequence = crossover(binary_sequence, round_info)
        encryption_key = get_pattern(identifier_key, key)[0]
        binary_sequence = bits_to_dna(
            group_bits(
                encryption_of_key(dna_to_bits(reverse_reshape(binary_sequence), utils.dna_base_to_two_bits_table), encryption_key)),
            utils.two_bits_to_dna_base_table)
        number_of_rounds -= 1

    return bin2str(dna_to_bits(binary_sequence, utils.dna_base_to_two_bits_table)).decode()


def main(encrypted_text):

    decrypted_file = open(decrypted_filename, "w")
    key_file = open(key_filename, "r")

    key = key_file.read()

    generate_pre_processing_tables()
    generate_mutation_tables()

    decrypted_text = dna_gdt(encrypted_text, key)

    decrypted_file.write(decrypted_text)

    decrypted_file.close()
    key_file.close()
    return decrypted_text
