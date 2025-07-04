import re

def clean_crtf(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    keywords = ['symsize', 'symlink', 'symthick', 'font', 'usetex']

    # Build one pattern that matches these keys
    # Handles cases like: , symsize=1 or symsize=1,
    pattern = re.compile(
        r'(,\s*)?(' + '|'.join(keywords) + r')\s*=\s*[^,]+(,)?'
    )

    cleaned_lines = []
    for line in lines:
        # Remove all matches iteratively
        cleaned = line
        while True:
            new_cleaned = pattern.sub(lambda m: ',' if m.group(1) and m.group(3) else '', cleaned)
            if new_cleaned == cleaned:
                break
            cleaned = new_cleaned
        # Clean up leftover extra commas
        cleaned = re.sub(r',\s*,', ',', cleaned)
        cleaned_lines.append(cleaned.strip() + '\n')

    with open(output_file, 'w') as f:
        f.writelines(cleaned_lines)

    print(f"Cleaned file written to: {output_file}")

# Example usage
clean_crtf('/orange/adamginsburg/w51/TaehwaYoo/regions/w51e_b3_taehwa_250703',
           '/orange/adamginsburg/w51/TaehwaYoo/regions/w51e_b3_taehwa_250703_cleaned.crtf')
clean_crtf('/orange/adamginsburg/w51/TaehwaYoo/regions/w51e_b6_taehwa_250703',
              '/orange/adamginsburg/w51/TaehwaYoo/regions/w51e_b6_taehwa_250703_cleaned.crtf')
clean_crtf('/orange/adamginsburg/w51/TaehwaYoo/regions/w51n_b3_taehwa_250703',
           '/orange/adamginsburg/w51/TaehwaYoo/regions/w51n_b3_taehwa_250703_cleaned.crtf')
clean_crtf('/orange/adamginsburg/w51/TaehwaYoo/regions/w51n_b6_taehwa_250703',
              '/orange/adamginsburg/w51/TaehwaYoo/regions/w51n_b6_taehwa_250703_cleaned.crtf')

