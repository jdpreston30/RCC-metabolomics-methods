#' Convert Metabolite Names to Sentence Case
#'
#' Applies IUPAC-compliant sentence case conventions to chemical compound names.
#' This function handles lowercase prefixes (locants and stereochemical descriptors)
#' and capitalizes the first letter of the actual chemical name following the prefix.
#'
#' @param name Character string of metabolite/chemical name
#' @param abbrev Character string indicating if name is an abbreviation ("Y" to skip conversion). Default NA.
#' @param preserve_markers Logical, whether to preserve superscript markers like *, ᶜ, ᴵ. Default TRUE.
#'
#' @details
#' According to IUPAC conventions, certain prefixes remain lowercase even at the
#' start of a sentence:
#' 
#' **Greek letter locants** (α, β, γ, δ, etc.):
#' - Remain lowercase
#' - Following letter is capitalized
#' - Example: "α-Ketoisocaproate" → "α-Ketoisocaproate"
#'
#' **Latin/Roman letter locants** (n-, p-, o-, m-, etc.):
#' - Written in italics (not handled in plain text)
#' - Remain lowercase
#' - Following letter is capitalized at sentence start
#' - Example: "n-Butyl iodide" → "n-Butyl iodide"
#'
#' **Stereochemical descriptors** (d-, l-, meso-, cis-, trans-, etc.):
#' - Written in italics (not handled in plain text)
#' - Remain lowercase
#' - Following letter is capitalized at sentence start
#' - Example: "d-Camphor" → "d-Camphor"
#'
#' **Positional/substituent descriptors** (N-, O-, S-, and numbered variants):
#' - Remain uppercase (these are element symbols)
#' - Following letter is capitalized at sentence start
#' - Example: "N-Acetyl-glucosaminate" → "N-acetyl-glucosaminate"
#'
#' **Abbreviations/acronyms**: Remain fully uppercase (e.g., GSSG, DHAP, ATP)
#'
#' @return Character string with proper sentence-case capitalization
#'
#' @references
#' International Union of Pure and Applied Chemistry (IUPAC) nomenclature conventions
#' Charlesworth Author Services: "Casing of Chemical Compounds 1: Rules for Lowercase Prefixes"
#'
#' @examples
#' metabolite_sentence_case("α-Ketoisocaproate")
#' # "α-Ketoisocaproate"
#'
#' metabolite_sentence_case("N-Acetyl-glucosaminate")
#' # "N-acetyl-glucosaminate"
#'
#' metabolite_sentence_case("1-O-16:0-LysoPC")
#' # "1-O-16:0-lysoPC"
#'
#' metabolite_sentence_case("β-methylphenylpyruvate")
#' # "β-Methylphenylpyruvate"
#'
#' @author Joshua D. Preston
#' @export
metabolite_sentence_case <- function(name, abbrev = NA, preserve_markers = TRUE) {
  
  if (is.na(name) || name == "") return(name)
  
  # If marked as abbreviation, return unchanged
  if (!is.na(abbrev) && abbrev == "Y") return(name)
  
  # Extract and preserve any trailing markers (*, ᶜ, ᴵ, etc.)
  marker <- ""
  if (preserve_markers) {
    marker_match <- regmatches(name, regexpr("[*ᶜᴵ]+$", name))
    if (length(marker_match) > 0) {
      marker <- marker_match
      name <- sub("[*ᶜᴵ]+$", "", name)
    }
  }
  
  # Greek letters that can appear as locants
  greek_letters <- c("α", "β", "γ", "δ", "ε", "ζ", "η", "θ", "κ", "λ", "μ", "ν", "ξ", "π", "ρ", "σ", "τ", "υ", "φ", "χ", "ψ", "ω")
  
  # Check if this contains consecutive capital letters (indicates abbreviation/acronym)
  # If so, preserve the original capitalization entirely
  if (grepl("[A-Z]{2,}", name)) {
    return(paste0(name, marker))
  }
  
  # Check if this is an all-caps abbreviation (2-6 characters, all uppercase, no hyphens after first char)
  if (grepl("^[A-Z]{2,6}$", name)) {
    return(paste0(name, marker))
  }
  
  # Pattern: starts with number-letter- (e.g., "1-O-", "2-OH-", "3-", "10,16-")
  # Capitalize first actual letter after the positional prefix
  if (grepl("^[0-9,]+-[A-Z]-", name)) {
    # Extract the prefix (e.g., "1-O-")
    prefix_match <- regmatches(name, regexpr("^[0-9,]+-[A-Z]-", name))
    prefix <- prefix_match[1]
    remainder <- substring(name, nchar(prefix) + 1)
    
    # Convert remainder to lowercase, then capitalize first letter
    remainder_lower <- tolower(remainder)
    remainder_cased <- paste0(toupper(substring(remainder_lower, 1, 1)), substring(remainder_lower, 2))
    
    return(paste0(prefix, remainder_cased, marker))
  }
  
  # Pattern: starts with number + Greek letter (e.g., "17α-estradiol")
  if (grepl("^[0-9,]+[αβγδεζηθκλμνξπρστυφχψω]-", name)) {
    # Extract number + Greek letter + hyphen
    prefix_match <- regmatches(name, regexpr("^[0-9,]+[αβγδεζηθκλμνξπρστυφχψω]-", name))
    prefix <- prefix_match[1]
    remainder <- substring(name, nchar(prefix) + 1)
    
    # Convert remainder to lowercase, then capitalize first letter
    remainder_lower <- tolower(remainder)
    remainder_cased <- paste0(toupper(substring(remainder_lower, 1, 1)), substring(remainder_lower, 2))
    
    return(paste0(prefix, remainder_cased, marker))
  }
  
  # Pattern: starts with number alone (e.g., "2-Aminobutanoate", "5-HTP", "3-(aminomethyl)indole")
  if (grepl("^[0-9,]+-", name)) {
    prefix_match <- regmatches(name, regexpr("^[0-9,]+-", name))
    prefix <- prefix_match[1]
    remainder <- substring(name, nchar(prefix) + 1)
    
    # Check if remainder starts with opening parenthesis
    if (grepl("^\\(", remainder)) {
      # Find the closing parenthesis and capitalize first letter after opening paren
      paren_content <- sub("^\\(([^)]+)\\)(.*)$", "\\1", remainder)
      after_paren <- sub("^\\([^)]+\\)(.*)$", "\\1", remainder)
      
      # Capitalize first letter of content inside parentheses
      paren_cased <- paste0(toupper(substring(paren_content, 1, 1)), tolower(substring(paren_content, 2)))
      
      # Lowercase and capitalize first letter after parentheses
      if (nchar(after_paren) > 0) {
        after_cased <- paste0(tolower(substring(after_paren, 1, 1)), tolower(substring(after_paren, 2)))
      } else {
        after_cased <- ""
      }
      
      return(paste0(prefix, "(", paren_cased, ")", after_cased, marker))
    } else {
      # Convert remainder to lowercase, then capitalize first letter
      remainder_lower <- tolower(remainder)
      remainder_cased <- paste0(toupper(substring(remainder_lower, 1, 1)), substring(remainder_lower, 2))
      
      return(paste0(prefix, remainder_cased, marker))
    }
  }
  
  # Pattern: starts with Greek letter followed by hyphen (e.g., "α-Keto...", "β-methyl...")
  first_char <- substring(name, 1, 1)
  if (first_char %in% greek_letters && substring(name, 2, 2) == "-") {
    # Greek letter stays lowercase, capitalize the letter after the hyphen
    remainder <- substring(name, 3)
    remainder_lower <- tolower(remainder)
    remainder_cased <- paste0(toupper(substring(remainder_lower, 1, 1)), substring(remainder_lower, 2))
    
    return(paste0(first_char, "-", remainder_cased, marker))
  }
  
  # Pattern: starts with capital letter descriptor (N-, O-, S-, L-, D- etc.)
  # These are element symbols or stereochemistry that should STAY UPPERCASE
  # Keep the descriptor capitalized and only process what follows
  if (grepl("^[A-Z]-", name)) {
    descriptor <- substring(name, 1, 2)  # E.g., "N-", "L-", "S-"
    remainder <- substring(name, 3)
    
    # Special handling: Process entire remainder but preserve capital letter descriptors
    # Look for patterns like -L-, -N-, -S-, -O-, -D- within the remainder and keep them uppercase
    result_parts <- c(descriptor)
    
    # Split by hyphens and process each part
    parts <- strsplit(remainder, "-")[[1]]
    first_word_found <- FALSE
    
    for (i in seq_along(parts)) {
      part <- parts[i]
      
      # Check if this part is a single capital letter (like "L", "N", "S", "O", "D")
      if (nchar(part) == 1 && grepl("^[A-Z]$", part)) {
        # Keep it uppercase and add hyphen
        result_parts <- c(result_parts, part)
        if (i < length(parts)) result_parts <- c(result_parts, "-")
      } else if (nchar(part) == 1 && part %in% greek_letters) {
        # Greek letter - keep lowercase
        result_parts <- c(result_parts, part)
        if (i < length(parts)) result_parts <- c(result_parts, "-")
      } else if (part == "(+)" || part == "(-)" || part == "(±)") {
        # Stereochemistry symbols - keep as-is
        result_parts <- c(result_parts, part)
        if (i < length(parts)) result_parts <- c(result_parts, "-")
      } else if (grepl("^\\(", part)) {
        # Starts with parenthesis - capitalize first letter after opening paren
        paren_match <- regmatches(part, regexpr("^\\([^)]*\\)", part))
        if (length(paren_match) > 0) {
          paren_content <- substring(paren_match, 2, nchar(paren_match) - 1)
          after_paren <- substring(part, nchar(paren_match) + 1)
          
          # Capitalize first letter inside parentheses
          paren_cased <- paste0(toupper(substring(paren_content, 1, 1)), tolower(substring(paren_content, 2)))
          
          # Handle text after closing paren - capitalize if first word
          if (nchar(after_paren) > 0 && !first_word_found) {
            after_cased <- paste0(toupper(substring(after_paren, 1, 1)), tolower(substring(after_paren, 2)))
            first_word_found <- TRUE
          } else if (nchar(after_paren) > 0) {
            after_cased <- tolower(after_paren)
          } else {
            after_cased <- ""
          }
          
          result_parts <- c(result_parts, paste0("(", paren_cased, ")", after_cased))
        } else {
          result_parts <- c(result_parts, part)
        }
        if (i < length(parts)) result_parts <- c(result_parts, "-")
      } else {
        # Regular part - capitalize only if it's the first actual word part
        if (!first_word_found) {
          # First word part - capitalize
          part_cased <- paste0(toupper(substring(part, 1, 1)), tolower(substring(part, 2)))
          first_word_found <- TRUE
        } else {
          # Subsequent parts - lowercase
          part_cased <- tolower(part)
        }
        result_parts <- c(result_parts, part_cased)
        if (i < length(parts)) result_parts <- c(result_parts, "-")
      }
    }
    
    return(paste0(paste(result_parts, collapse = ""), marker))
  }
  
  # Pattern: starts with lowercase descriptor (d-, l-, meso-, cis-, trans-, etc.)
  # These stay lowercase, capitalize what follows
  stereo_descriptors <- c("d-", "l-", "cis-", "trans-", "meso-", "rac-", "n-", "o-", "p-", "m-", "sec-", "tert-")
  for (descriptor in stereo_descriptors) {
    if (grepl(paste0("^", descriptor), name, ignore.case = TRUE)) {
      # Only apply if it's actually lowercase in the input
      if (substring(name, 1, nchar(descriptor)) == tolower(substring(name, 1, nchar(descriptor)))) {
        desc_len <- nchar(descriptor)
        remainder <- substring(name, desc_len + 1)
        remainder_lower <- tolower(remainder)
        remainder_cased <- paste0(toupper(substring(remainder_lower, 1, 1)), substring(remainder_lower, 2))
        
        return(paste0(tolower(descriptor), remainder_cased, marker))
      }
    }
  }
  
  # Default: standard sentence case (capitalize first letter, rest lowercase)
  name_lower <- tolower(name)
  result <- paste0(toupper(substring(name_lower, 1, 1)), substring(name_lower, 2))
  
  return(paste0(result, marker))
}
