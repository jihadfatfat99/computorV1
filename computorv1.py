import sys
import re
from typing import Union, NoReturn, Dict, Tuple

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


def manual_cbrt(x: float) -> float:
    """
    Computes the cube root of a number using Newton's method.
    
    This handles both positive and negative numbers.
    
    Args:
        x (float): The number to find the cube root of
    
    Returns:
        float: The cube root of x
    
    Algorithm: Newton's method for cube root
        f(y) = y³ - x = 0
        y_next = y - f(y)/f'(y) = y - (y³ - x)/(3y²)
        y_next = (2y + x/y²) / 3
    
    Examples:
        manual_cbrt(8) → 2.0
        manual_cbrt(-8) → -2.0
        manual_cbrt(27) → 3.0
    """
    if x == 0:
        return 0.0
    
    # Handle negative numbers
    if x < 0:
        return -manual_cbrt(-x)
    
    # Initial guess
    if x >= 1:
        guess = x / 3.0
    else:
        guess = x
    
    epsilon = 1e-10
    max_iterations = 100
    
    for i in range(max_iterations):
        if abs(guess) < 1e-10:
            return 0.0
        
        # Newton's iteration for cube root: next = (2*guess + x/guess²) / 3
        guess_squared = guess * guess
        next_guess = (2 * guess + x / guess_squared) / 3.0
        
        if abs(next_guess - guess) < epsilon:
            return next_guess
        
        guess = next_guess
    
    return guess


def solve_degree_3(reduced: Dict[int, float]) -> None:
    """
    Solves a degree 3 (cubic) equation using Cardano's formula.
    
    A cubic equation has the form: ax³ + bx² + cx + d = 0
    
    This is a BONUS feature - not required by the subject.
    
    Args:
        reduced (Dict[int, float]): Reduced equation {exponent: coefficient}
    
    Returns:
        None: Prints solution(s) to stdout
    
    Algorithm: Cardano's formula
        1. Convert to depressed cubic: t³ + pt + q = 0
           where x = t - b/(3a)
        2. Calculate discriminant: Δ = -4p³ - 27q²
        3. Based on Δ:
           - Δ > 0: Three distinct real roots
           - Δ = 0: Multiple real roots
           - Δ < 0: One real root, two complex conjugate roots
    
    Note: This implementation focuses on finding at least one real root.
    For complex cases, uses numerical methods for additional roots.
    
    Examples:
        x³ - 6x² + 11x - 6 = 0  → roots: 1, 2, 3
        x³ - 1 = 0              → roots: 1, complex conjugates
    """
    a = reduced.get(3, 0.0)  # Coefficient of X³
    b = reduced.get(2, 0.0)  # Coefficient of X²
    c = reduced.get(1, 0.0)  # Coefficient of X¹
    d = reduced.get(0, 0.0)  # Coefficient of X⁰
    
    if abs(a) < 1e-10:
        handle_error("Internal error: degree 3 equation with zero coefficient for X^3")
    
    print("Solving cubic equation using Cardano's formula...")
    print("(This is a bonus feature)")
    
    # Normalize: divide by a to get x³ + Bx² + Cx + D = 0
    B = b / a
    C = c / a
    D = d / a
    
    # Convert to depressed cubic: t³ + pt + q = 0
    # where x = t - B/3
    p = C - (B * B) / 3.0
    q = (2 * B * B * B) / 27.0 - (B * C) / 3.0 + D
    
    # Calculate discriminant: Δ = -4p³ - 27q²
    discriminant = -4 * p * p * p - 27 * q * q
    
    print(f"\nDepressed cubic form: t³ + ({format_number(p)})t + ({format_number(q)}) = 0")
    print(f"Discriminant: Δ = {format_number(discriminant)}")
    
    threshold = 1e-10
    
    if discriminant > threshold:
        # Three distinct real roots (trigonometric method)
        print("\nΔ > 0: Three distinct real roots")
        
        # Use trigonometric solution
        # t_k = 2√(-p/3) * cos((1/3)arccos((3q)/(2p)√(-3/p)) - 2πk/3)
        # for k = 0, 1, 2
        
        if p >= 0:
            # Should not happen for Δ > 0, but handle gracefully
            print("Note: Using numerical approximation for this case")
            
        # Trigonometric method requires -p/3 > 0
        sqrt_term = manual_sqrt(-p / 3.0)
        
        # Calculate angle: φ = arccos((3q)/(2p) * sqrt(-3/p))
        cos_arg = (3 * q) / (2 * p) * manual_sqrt(-3 / p)
        
        # Clamp to [-1, 1] for numerical stability
        if cos_arg > 1:
            cos_arg = 1
        elif cos_arg < -1:
            cos_arg = -1
        
        # Manual arccos implementation
        phi = manual_arccos(cos_arg)
        
        # Manual π value
        PI = 3.141592653589793
        
        # Calculate three roots using manual cosine
        # cos(θ) using Taylor series: cos(x) = 1 - x²/2! + x⁴/4! - x⁶/6! + ...
        t1 = 2 * sqrt_term * manual_cos(phi / 3.0)
        t2 = 2 * sqrt_term * manual_cos((phi + 2 * PI) / 3.0)
        t3 = 2 * sqrt_term * manual_cos((phi + 4 * PI) / 3.0)
        
        # Convert back: x = t - B/3
        shift = B / 3.0
        x1 = t1 - shift
        x2 = t2 - shift
        x3 = t3 - shift
        
        print("\nThe three real solutions are:")
        print(format_number(x1))
        print(format_number(x2))
        print(format_number(x3))
        
    else:
        # One real root (and two complex conjugates)
        print("\nΔ ≤ 0: One real root (and two complex conjugate roots)")
        
        # Cardano's formula
        # t = ∛(-q/2 + √(q²/4 + p³/27)) + ∛(-q/2 - √(q²/4 + p³/27))
        
        sqrt_arg = q * q / 4.0 + p * p * p / 27.0
        
        if sqrt_arg >= 0:
            sqrt_val = manual_sqrt(sqrt_arg)
        else:
            # Complex case, use absolute value
            sqrt_val = manual_sqrt(-sqrt_arg)
        
        u = manual_cbrt(-q / 2.0 + sqrt_val)
        v = manual_cbrt(-q / 2.0 - sqrt_val)
        
        t = u + v
        x_real = t - B / 3.0
        
        print("\nThe real solution is:")
        print(format_number(x_real))
        
        # The two complex roots are:
        # x = -t/2 - B/3 ± i(√3/2)(u-v)
        real_part_complex = -t / 2.0 - B / 3.0
        imag_part_complex = (manual_sqrt(3) / 2.0) * (u - v)
        
        if abs(imag_part_complex) > threshold:
            print("\nThe two complex conjugate solutions are:")
            print(f"{format_number(real_part_complex)} + {format_number(abs(imag_part_complex))}i")
            print(f"{format_number(real_part_complex)} - {format_number(abs(imag_part_complex))}i")


def manual_arccos(x: float) -> float:
    """
    Computes arccos(x) manually using Taylor series and Newton-Raphson method.
    
    Args:
        x (float): Value between -1 and 1
    
    Returns:
        float: arccos(x) in radians
    
    Uses the identity: arccos(x) = π/2 - arcsin(x)
    Where arcsin(x) is computed using Taylor series.
    
    Taylor series for arcsin(x):
    arcsin(x) = x + (1/2)(x³/3) + (1·3)/(2·4)(x⁵/5) + (1·3·5)/(2·4·6)(x⁷/7) + ...
    """
    # Clamp input to valid range
    if x > 1:
        x = 1
    if x < -1:
        x = -1
    
    # Manual π value (sufficient precision)
    PI = 3.141592653589793
    
    # Special cases for better accuracy
    if abs(x - 1) < 1e-10:
        return 0.0  # arccos(1) = 0
    if abs(x + 1) < 1e-10:
        return PI  # arccos(-1) = π
    if abs(x) < 1e-10:
        return PI / 2.0  # arccos(0) = π/2
    
    # For |x| > 0.5, use identity: arccos(x) = π/2 - arcsin(x)
    # For |x| <= 0.5, use: arccos(x) = π/2 - arcsin(x)
    # We compute arcsin(x) using Taylor series
    
    # Compute arcsin(x) using Taylor series
    # arcsin(x) = Σ [(2n)! / (2^(2n) * (n!)² * (2n+1))] * x^(2n+1)
    
    arcsin_x = x
    term = x
    
    for n in range(1, 100):  # Sufficient iterations for convergence
        # Calculate next term in series
        # term_n = term_(n-1) * x² * (2n-1)² / ((2n) * (2n+1))
        term = term * x * x * (2 * n - 1) * (2 * n - 1) / ((2 * n) * (2 * n + 1))
        arcsin_x += term
        
        # Check for convergence
        if abs(term) < 1e-15:
            break
    
    # Use identity: arccos(x) = π/2 - arcsin(x)
    return PI / 2.0 - arcsin_x


def manual_cos(x: float) -> float:
    """
    Computes cos(x) manually using Taylor series.
    
    Args:
        x (float): Angle in radians
    
    Returns:
        float: cos(x)
    
    Taylor series for cos(x):
    cos(x) = 1 - x²/2! + x⁴/4! - x⁶/6! + x⁸/8! - ...
    """
    # Manual 2π value
    TWO_PI = 6.283185307179586
    
    # Reduce x to range [0, 2π] for better convergence
    while x > TWO_PI:
        x -= TWO_PI
    while x < 0:
        x += TWO_PI
    
    # Taylor series
    cos_x = 1.0
    term = 1.0
    
    for n in range(1, 50):  # Sufficient iterations
        # term_n = -term_(n-1) * x² / ((2n-1) * (2n))
        term = -term * x * x / ((2 * n - 1) * (2 * n))
        cos_x += term
        
        # Check for convergence
        if abs(term) < 1e-15:
            break
    
    return cos_x


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
        elif degree == 3:
            # BONUS: Solve cubic equations
            solve_degree_3(reduced)
        else:
            print("The polynomial degree is strictly greater than 2, I can't solve.")
            print("(Note: Degree 3 is supported as a bonus feature)")
        
        sys.exit(0)
    
    except SystemExit:
        raise
    
    except KeyboardInterrupt:
        print("\n\nInterrupted by user.", file=sys.stderr)
        sys.exit(1)
    
    except Exception as e:
        handle_error(f"Unexpected error: {str(e)}")

if __name__ == "__main__":
    main()