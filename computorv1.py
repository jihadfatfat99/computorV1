#!/usr/bin/env python3
"""
Computor v1 - Polynomial Equation Solver
42 School Project - Cryptography and Math Specialization

This program solves polynomial equations of degree 2 or lower.
It displays the reduced form, polynomial degree, and solution(s).

Usage:
    ./computorv1 "equation"
    ./computorv1  (reads from STDIN)

Examples:
    ./computorv1 "5 * X^0 + 4 * X^1 - 9.3 * X^2 = 1 * X^0"
    echo "X^2 + 5X + 3 = 0" | ./computorv1
"""

import sys
import re
from typing import Union, NoReturn, Dict, Tuple


# ============================================================================
# PART 1: UTILITY AND ERROR FUNCTIONS
# ============================================================================

def handle_error(message: str) -> NoReturn:
    """Prints an error message to stderr and exits the program."""
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(1)


def format_number(n: Union[int, float, None]) -> str:
    """Formats a number for clean display in output."""
    if n is None:
        return "0"
    
    if isinstance(n, float):
        if n != n:  # NaN check
            return "nan"
        if n == float('inf'):
            return "inf"
        if n == float('-inf'):
            return "-inf"
    
    num = float(n)
    
    if num == 0:
        return "0"
    
    if abs(num) < 1e-10:
        return "0"
    
    num = round(num, 6)
    
    if num == int(num):
        return str(int(num))
    
    formatted = f"{num:.6f}".rstrip('0').rstrip('.')
    return formatted


def gcd(a: Union[int, float], b: Union[int, float]) -> int:
    """Computes the Greatest Common Divisor using Euclidean algorithm."""
    a = abs(int(a))
    b = abs(int(b))
    
    if a == 0 and b == 0:
        return 1
    
    while b != 0:
        a, b = b, a % b
    
    return a


def simplify_fraction(numerator: Union[int, float], 
                     denominator: Union[int, float]) -> str:
    """Reduces a fraction to its irreducible form."""
    if denominator == 0:
        handle_error("Division by zero in fraction simplification")
    
    if abs(numerator) < 1e-10:
        return "0"
    
    sign = 1
    if (numerator < 0) != (denominator < 0):
        sign = -1
    
    num = abs(numerator)
    den = abs(denominator)
    
    precision = 1000000
    num_int = round(num * precision)
    den_int = round(den * precision)
    
    divisor = gcd(num_int, den_int)
    
    if divisor > 0:
        num_int //= divisor
        den_int //= divisor
    
    num_int *= sign
    
    if den_int == 1:
        return str(num_int)
    
    if den_int > 10000:
        result = numerator / denominator
        return format_number(result)
    
    return f"{num_int}/{den_int}"


# ============================================================================
# PART 2: PARSING UTILITIES
# ============================================================================

def validate_equation(equation_str: str) -> bool:
    """Validates the equation string for correct syntax."""
    if not equation_str or equation_str.strip() == "":
        handle_error("Empty equation provided")
    
    equals_count = equation_str.count('=')
    if equals_count == 0:
        handle_error("Missing '=' sign in equation")
    if equals_count > 1:
        handle_error("Multiple '=' signs found in equation")
    
    parts = equation_str.split('=')
    left_side = parts[0].strip()
    right_side = parts[1].strip()
    
    if not left_side:
        handle_error("Left side of equation is empty")
    if not right_side:
        handle_error("Right side of equation is empty")
    
    allowed_pattern = re.compile(r'^[0-9X\^\*\+\-=\s\.]+$')
    
    if not allowed_pattern.match(equation_str):
        invalid_chars = set()
        for char in equation_str:
            if not re.match(r'[0-9X\^\*\+\-=\s\.]', char):
                invalid_chars.add(char)
        
        if invalid_chars:
            handle_error(f"Invalid characters found: {', '.join(sorted(invalid_chars))}")
        else:
            handle_error("Invalid equation format")
    
    if re.search(r'[\+\-\*\^]{3,}', equation_str):
        handle_error("Invalid operator sequence")
    
    if re.match(r'^[\+\*\^=]', equation_str.strip()):
        handle_error("Equation cannot start with an operator")
    if re.search(r'[\+\-\*\^]$', equation_str.strip()):
        handle_error("Equation cannot end with an operator")
    
    return True


def parse_freeform_term(term: str) -> Tuple[float, int]:
    """Parses a single term and extracts coefficient and exponent."""
    term = term.strip()
    
    if not term:
        return (0.0, 0)
    
    coefficient = 1.0
    exponent = 0
    
    # Pattern 1: Standard format "a * X^p"
    standard_pattern = r'^([+-]?\d*\.?\d+)\s*\*\s*X\s*\^\s*(\d+)$'
    match = re.match(standard_pattern, term, re.IGNORECASE)
    if match:
        coefficient = float(match.group(1))
        exponent = int(match.group(2))
        return (coefficient, exponent)
    
    # Pattern 2: Free form "aX^p"
    freeform_pattern = r'^([+-]?\d*\.?\d+)\s*\*?\s*X\s*\^\s*(\d+)$'
    match = re.match(freeform_pattern, term, re.IGNORECASE)
    if match:
        coef_str = match.group(1)
        coefficient = float(coef_str) if coef_str else 1.0
        exponent = int(match.group(2))
        return (coefficient, exponent)
    
    # Pattern 3: Just "X^p"
    x_power_pattern = r'^([+-]?)X\s*\^\s*(\d+)$'
    match = re.match(x_power_pattern, term, re.IGNORECASE)
    if match:
        sign = match.group(1)
        coefficient = -1.0 if sign == '-' else 1.0
        exponent = int(match.group(2))
        return (coefficient, exponent)
    
    # Pattern 4: Free form "aX"
    freeform_x_pattern = r'^([+-]?\d*\.?\d+)\s*\*?\s*X$'
    match = re.match(freeform_x_pattern, term, re.IGNORECASE)
    if match:
        coef_str = match.group(1)
        coefficient = float(coef_str) if coef_str else 1.0
        exponent = 1
        return (coefficient, exponent)
    
    # Pattern 5: Just "X"
    if re.match(r'^[+-]?X$', term, re.IGNORECASE):
        coefficient = -1.0 if term.strip().startswith('-') else 1.0
        exponent = 1
        return (coefficient, exponent)
    
    # Pattern 6: Just a number
    number_pattern = r'^([+-]?\d+\.?\d*)$'
    match = re.match(number_pattern, term)
    if match:
        coefficient = float(match.group(1))
        exponent = 0
        return (coefficient, exponent)
    
    handle_error(f"Invalid term format: '{term}'")


def extract_terms(side_str: str) -> Dict[int, float]:
    """Extracts all terms from one side of the equation."""
    terms_dict: Dict[int, float] = {}
    
    side_str = side_str.strip()
    
    if not side_str:
        return terms_dict
    
    side_str = re.sub(r'\s*\+\s*', ' +', side_str)
    side_str = re.sub(r'\s*-\s*', ' -', side_str)
    
    parts = side_str.split()
    
    terms = []
    current_term = ""
    
    for part in parts:
        if part in ['+', '-']:
            if current_term:
                terms.append(current_term)
            current_term = part
        elif part.startswith('+') or part.startswith('-'):
            if current_term:
                terms.append(current_term)
            current_term = part
        else:
            if current_term and not current_term[-1].isspace():
                current_term += " " + part
            else:
                current_term += part
    
    if current_term:
        terms.append(current_term)
    
    if not terms:
        terms = [side_str]
    
    for term in terms:
        term = term.strip()
        if not term or term in ['+', '-']:
            continue
        
        try:
            coefficient, exponent = parse_freeform_term(term)
            
            if exponent in terms_dict:
                terms_dict[exponent] += coefficient
            else:
                terms_dict[exponent] = coefficient
        
        except Exception as e:
            handle_error(f"Failed to parse term '{term}': {str(e)}")
    
    return terms_dict


def parse_input(equation_str: str) -> Tuple[Dict[int, float], Dict[int, float]]:
    """Parses the complete equation into left and right terms."""
    validate_equation(equation_str)
    
    parts = equation_str.split('=')
    left_side = parts[0].strip()
    right_side = parts[1].strip()
    
    left_terms = extract_terms(left_side)
    right_terms = extract_terms(right_side)
    
    return (left_terms, right_terms)


# ============================================================================
# PART 3: EQUATION REDUCTION AND DEGREE
# ============================================================================

def reduce_equation(left_terms: Dict[int, float], 
                   right_terms: Dict[int, float]) -> Dict[int, float]:
    """Reduces the equation by moving all terms to the left side."""
    reduced: Dict[int, float] = {}
    
    all_exponents = set(left_terms.keys()) | set(right_terms.keys())
    
    for exp in all_exponents:
        left_coef = left_terms.get(exp, 0.0)
        right_coef = right_terms.get(exp, 0.0)
        reduced[exp] = left_coef - right_coef
    
    if not reduced:
        reduced[0] = 0.0
    
    return reduced


def get_polynomial_degree(reduced: Dict[int, float]) -> int:
    """Determines the degree of the polynomial equation."""
    if not reduced:
        return 0
    
    threshold = 1e-10
    
    max_degree = 0
    has_nonzero = False
    
    for exp, coef in reduced.items():
        if abs(coef) > threshold:
            has_nonzero = True
            if exp > max_degree:
                max_degree = exp
    
    if not has_nonzero:
        return 0
    
    return max_degree


# ============================================================================
# PART 4: OUTPUT DISPLAY
# ============================================================================

def print_reduced_form(reduced: Dict[int, float]) -> None:
    """Prints the reduced form of the equation in the required format."""
    print("Reduced form: ", end="")
    
    if not reduced:
        print("0 * X^0 = 0")
        return
    
    threshold = 1e-10
    all_zero = all(abs(coef) < threshold for coef in reduced.values())
    
    if all_zero:
        print("0 * X^0 = 0")
        return
    
    sorted_exponents = sorted(reduced.keys())
    
    terms = []
    
    for exp in sorted_exponents:
        coef = reduced[exp]
        coef_str = format_number(coef)
        term = f"{coef_str} * X^{exp}"
        terms.append((coef, term))
    
    if terms:
        first_coef, first_term = terms[0]
        print(first_term, end="")
        
        for i in range(1, len(terms)):
            coef, term = terms[i]
            
            if coef >= 0:
                print(f" + {term}", end="")
            else:
                print(f" - {term[1:]}", end="")
        
        print(" = 0")


def print_steps(reduced: Dict[int, float], degree: int) -> None:
    """Prints intermediate calculation steps for solving the equation."""
    print("\n--- Solution Steps ---")
    
    if degree == 0:
        c = reduced.get(0, 0.0)
        print(f"This is a constant equation: {format_number(c)} = 0")
        
        if abs(c) < 1e-10:
            print("Since 0 = 0 is always true, any real number is a solution.")
        else:
            print(f"Since {format_number(c)} ≠ 0, there is no solution.")
    
    elif degree == 1:
        a = reduced.get(1, 0.0)
        b = reduced.get(0, 0.0)
        
        print(f"This is a linear equation: {format_number(a)}X + {format_number(b)} = 0")
        print(f"Using formula: x = -b/a")
        print(f"Where a = {format_number(a)}, b = {format_number(b)}")
        print(f"Calculating: x = -({format_number(b)}) / {format_number(a)}")
        print(f"           x = {format_number(-b)} / {format_number(a)}")
        
        solution = -b / a
        print(f"           x = {format_number(solution)}")
    
    elif degree == 2:
        a = reduced.get(2, 0.0)
        b = reduced.get(1, 0.0)
        c = reduced.get(0, 0.0)
        
        print(f"This is a quadratic equation: {format_number(a)}X² + {format_number(b)}X + {format_number(c)} = 0")
        print(f"Coefficients: a = {format_number(a)}, b = {format_number(b)}, c = {format_number(c)}")
        print(f"\nUsing quadratic formula: x = (-b ± √Δ) / 2a")
        print(f"Where Δ = b² - 4ac (discriminant)")
        
        discriminant = b * b - 4 * a * c
        
        print(f"\nCalculating discriminant:")
        print(f"  Δ = ({format_number(b)})² - 4({format_number(a)})({format_number(c)})")
        print(f"  Δ = {format_number(b * b)} - {format_number(4 * a * c)}")
        print(f"  Δ = {format_number(discriminant)}")
        
        if discriminant > 0:
            print(f"\nΔ > 0: Two distinct real solutions exist")
        elif abs(discriminant) < 1e-10:
            print(f"\nΔ = 0: One repeated real solution exists")
        else:
            print(f"\nΔ < 0: Two complex conjugate solutions exist")
    
    else:
        print(f"This is a degree {degree} polynomial equation.")
        print("Equations of degree higher than 2 cannot be solved by this program.")
    
    print("--- End Steps ---\n")


# ============================================================================
# PART 5: MATHEMATICAL SOLVERS
# ============================================================================

def manual_sqrt(x: float) -> float:
    """Computes the square root using Newton's method (Babylonian method)."""
    if x < 0:
        handle_error(f"Cannot compute square root of negative number: {x}")
    
    if x == 0 or abs(x) < 1e-10:
        return 0.0
    
    if abs(x - 1.0) < 1e-10:
        return 1.0
    
    if x >= 1:
        guess = x / 2.0
    else:
        guess = x
    
    epsilon = 1e-10
    max_iterations = 100
    
    for i in range(max_iterations):
        if guess == 0:
            return 0.0
        
        next_guess = (guess + x / guess) / 2.0
        
        if abs(next_guess - guess) < epsilon:
            return next_guess
        
        guess = next_guess
    
    return guess


def solve_degree_0(reduced: Dict[int, float]) -> None:
    """Solves a degree 0 (constant) equation."""
    c = reduced.get(0, 0.0)
    
    threshold = 1e-10
    
    if abs(c) < threshold:
        print("Any real number is a solution.")
    else:
        print("No solution.")


def solve_degree_1(reduced: Dict[int, float]) -> None:
    """Solves a degree 1 (linear) equation."""
    a = reduced.get(1, 0.0)
    b = reduced.get(0, 0.0)
    
    if abs(a) < 1e-10:
        handle_error("Internal error: degree 1 equation with zero coefficient for X^1")
    
    solution = -b / a
    
    print("The solution is:")
    
    fraction_str = simplify_fraction(-b, a)
    
    if '/' in fraction_str and len(fraction_str) <= 10:
        print(f"{fraction_str}")
        if abs(solution - int(solution)) > 1e-10:
            print(f"({format_number(solution)})")
    else:
        print(format_number(solution))


def solve_degree_2(reduced: Dict[int, float]) -> None:
    """Solves a degree 2 (quadratic) equation."""
    a = reduced.get(2, 0.0)
    b = reduced.get(1, 0.0)
    c = reduced.get(0, 0.0)
    
    if abs(a) < 1e-10:
        handle_error("Internal error: degree 2 equation with zero coefficient for X^2")
    
    discriminant = b * b - 4 * a * c
    
    threshold = 1e-10
    
    if discriminant > threshold:
        print("Discriminant is strictly positive, the two solutions are:")
        
        sqrt_discriminant = manual_sqrt(discriminant)
        
        solution1 = (-b + sqrt_discriminant) / (2 * a)
        solution2 = (-b - sqrt_discriminant) / (2 * a)
        
        print(format_number(solution1))
        print(format_number(solution2))
    
    elif abs(discriminant) <= threshold:
        print("Discriminant is zero, the solution is:")
        
        solution = -b / (2 * a)
        print(format_number(solution))
    
    else:
        print("Discriminant is strictly negative, the two complex solutions are:")
        
        real_part = -b / (2 * a)
        
        abs_discriminant = -discriminant
        sqrt_abs_discriminant = manual_sqrt(abs_discriminant)
        imag_part = sqrt_abs_discriminant / (2 * a)
        
        real_frac = simplify_fraction(-b, 2 * a)
        imag_frac = simplify_fraction(sqrt_abs_discriminant, 2 * a)
        
        if '/' in real_frac and '/' in imag_frac:
            print(f"{real_frac} + {imag_frac}i")
            print(f"{real_frac} - {imag_frac}i")
        else:
            real_str = format_number(real_part)
            imag_str = format_number(abs(imag_part))
            
            if imag_part >= 0:
                print(f"{real_str} + {imag_str}i")
                print(f"{real_str} - {imag_str}i")
            else:
                print(f"{real_str} - {imag_str}i")
                print(f"{real_str} + {imag_str}i")


# ============================================================================
# PART 6: MAIN PROGRAM
# ============================================================================

def main() -> None:
    """Main program orchestrator - coordinates the entire equation solving process."""
    try:
        # Read input
        if len(sys.argv) > 1:
            equation_str = sys.argv[1]
        else:
            try:
                equation_str = input()
            except EOFError:
                handle_error("No equation provided")
        
        if not equation_str or equation_str.strip() == "":
            handle_error("Empty equation provided")
        
        # Validate equation
        validate_equation(equation_str)
        
        # Parse equation
        left_terms, right_terms = parse_input(equation_str)
        
        # Reduce equation
        reduced = reduce_equation(left_terms, right_terms)
        
        # Display reduced form
        print_reduced_form(reduced)
        
        # Determine polynomial degree
        degree = get_polynomial_degree(reduced)
        
        # Display polynomial degree
        print(f"Polynomial degree: {degree}")
        
        # Solve based on degree
        if degree == 0:
            solve_degree_0(reduced)
        elif degree == 1:
            solve_degree_1(reduced)
        elif degree == 2:
            solve_degree_2(reduced)
        else:
            print("The polynomial degree is strictly greater than 2, I can't solve.")
        
        sys.exit(0)
    
    except SystemExit:
        raise
    
    except KeyboardInterrupt:
        print("\n\nInterrupted by user.", file=sys.stderr)
        sys.exit(1)
    
    except Exception as e:
        handle_error(f"Unexpected error: {str(e)}")


# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    main()